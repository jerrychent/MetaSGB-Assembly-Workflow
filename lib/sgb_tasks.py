import ast
import os
import shutil
from pathlib import Path

import pandas as pd


"""Task builders for the MetaSGB assembly workflow."""


def _bool_flag(enabled, flag):
    """Return a CLI flag only when the option is enabled."""
    return flag if enabled else ""


def run_kneaddata(
    workflow,
    input_files,
    output_dir,
    db_path,
    suffix_r1="_1.fastq.gz",
    suffix_r2="_2.fastq.gz",
    threads=4,
    scheduler_cores=None,
    mem_mb=10240,
    partition=None,
    sequencer_source="TruSeq3",
    run_fastqc_start=True,
    run_fastqc_end=True,
    extra_args="",
):
    """Create one KneadData task per paired-end sample."""
    scheduler_cores = scheduler_cores or threads
    clean_files_info = []
    fastqc_flags = " ".join(
        flag
        for flag in [
            _bool_flag(run_fastqc_start, "--run-fastqc-start"),
            _bool_flag(run_fastqc_end, "--run-fastqc-end"),
        ]
        if flag
    )

    for r1 in input_files:
        dirname = os.path.dirname(r1)
        basename = os.path.basename(r1)
        if not basename.endswith(suffix_r1):
            continue

        sample = basename.rsplit(suffix_r1, 1)[0]
        r2 = os.path.join(dirname, basename.replace(suffix_r1, suffix_r2))
        if not os.path.exists(r2):
            raise FileNotFoundError(f"R2 file not found: {r2}")

        sample_out_dir = os.path.join(output_dir, sample)
        clean_r1 = os.path.join(sample_out_dir, f"{sample}_paired_1.fastq")
        clean_r2 = os.path.join(sample_out_dir, f"{sample}_paired_2.fastq")
        log_file = os.path.join(sample_out_dir, f"{sample}_kneaddata.log")

        # The active environment supplies KneadData; the database is selected per run.
        cmd = (
            f"/bin/bash -c '"
            f"set -euo pipefail; "
            f"mkdir -p [args[0]]; "
            f"exec > [args[7]] 2>&1; "
            f"kneaddata --input1 [depends[0]] --input2 [depends[1]] "
            f"--output [args[0]] --output-prefix [args[1]] "
            f"-t [args[2]] -db [args[3]] --bowtie2-options=\"-p [args[2]]\" "
            f"--sequencer-source [args[4]] [args[5]] [args[6]] "
            f"'"
        )

        workflow.add_task_gridable(
            cmd,
            depends=[r1, r2],
            targets=[clean_r1, clean_r2],
            args=[sample_out_dir, sample, threads, db_path, sequencer_source, fastqc_flags, extra_args, log_file],
            cores=scheduler_cores,
            mem=mem_mb,
            partition=partition,
            time=2880,
            name=f"knead_{sample}",
        )
        clean_files_info.append((sample, clean_r1, clean_r2))

    return clean_files_info


def run_assembly(workflow, clean_files_info, output_base_dir, threads=16, scheduler_cores=None, mem_mb=200000, partition=None, extra_args=""):
    """Create one metaSPAdes assembly task per sample."""
    scheduler_cores = scheduler_cores or threads
    assembly_info = []

    for sample, r1, r2 in clean_files_info:
        sample_out_dir = os.path.join(output_base_dir, sample)
        contigs = os.path.join(sample_out_dir, "contigs.fasta")
        log_file = os.path.join(sample_out_dir, "spades_console.log")

        cmd = (
            f"/bin/bash -c '"
            f"set -euo pipefail; "
            f"mkdir -p [args[0]]; "
            f"exec > [args[4]] 2>&1; "
            f"metaspades.py -1 [depends[0]] -2 [depends[1]] "
            f"-o [args[0]] -t [args[1]] --memory [args[2]] [args[3]] "
            f"'"
        )

        workflow.add_task_gridable(
            cmd,
            depends=[r1, r2],
            targets=contigs,
            args=[sample_out_dir, threads, int(mem_mb / 1024), extra_args, log_file],
            cores=scheduler_cores,
            mem=mem_mb,
            partition=partition,
            time=10080,
            name=f"assembly_{sample}",
        )
        assembly_info.append((sample, contigs))

    return assembly_info


def run_binning_workflow(
    workflow,
    assembly_info,
    clean_files_map,
    output_base_dir,
    threads=16,
    scheduler_cores=None,
    mem_mb=32000,
    partition=None,
    min_contig=1500,
    bowtie2_extra_args="--very-sensitive-local",
    metabat_extra_args="",
):
    """Create Bowtie2, MetaBAT2, and CheckM tasks for each assembly."""
    scheduler_cores = scheduler_cores or threads
    all_checkm_stats = []
    dir_binning = os.path.join(output_base_dir, "04_Binning")
    dir_checkm = os.path.join(output_base_dir, "05_Checkm")

    for sample, contigs in assembly_info:
        if sample not in clean_files_map:
            print(f"Warning: Clean data for {sample} not found, skipping binning.")
            continue
        r1, r2 = clean_files_map[sample]

        assembly_dir = os.path.dirname(contigs)
        index_base = os.path.join(assembly_dir, "contigs_index")
        index_done = index_base + ".done"
        sorted_bam = os.path.join(assembly_dir, f"{sample}.sorted.bam")
        sorted_bam_index = sorted_bam + ".bai"
        log_build = os.path.join(assembly_dir, "bowtie2_build.log")

        # Track Bowtie2 index completion with a sentinel to support .bt2 and .bt2l indexes.
        workflow.add_task_gridable(
            f"/bin/bash -c '"
            f"set -euo pipefail; "
            f"mkdir -p [args[2]]; "
            f"bowtie2-build --threads [args[0]] [depends[0]] [args[1]] > [args[3]] 2>&1; "
            f"touch [targets[0]]"
            f"'",
            depends=contigs,
            targets=index_done,
            args=[threads, index_base, assembly_dir, log_build],
            cores=scheduler_cores,
            mem=mem_mb,
            partition=partition,
            time=1440,
            name=f"index_{sample}",
        )

        # Local very-sensitive mapping improves placement for fragmented MAG contigs.
        log_map = os.path.join(assembly_dir, "bowtie2_mapping.log")
        map_mem = max(mem_mb, 32000)
        workflow.add_task_gridable(
            f"/bin/bash -c '"
            f"set -euo pipefail; "
            f"bowtie2 -p [args[0]] [args[2]] -x [args[1]] -1 [depends[0]] -2 [depends[1]] | "
            f"samtools sort -@ [args[0]] -o [targets[0]] -; "
            f"samtools index [targets[0]]"
            f"' > [args[3]] 2>&1",
            depends=[r1, r2, index_done],
            targets=[sorted_bam, sorted_bam_index],
            args=[threads, index_base, bowtie2_extra_args, log_map],
            cores=scheduler_cores,
            mem=map_mem,
            partition=partition,
            time=2880,
            name=f"map_{sample}",
        )

        sample_bin_dir = os.path.join(dir_binning, sample)
        depth_file = os.path.join(sample_bin_dir, f"{sample}.depth.txt")
        metabat_done = os.path.join(sample_bin_dir, "metabat.done")
        log_metabat = os.path.join(sample_bin_dir, "metabat.log")

        workflow.add_task_gridable(
            f"/bin/bash -c '"
            f"set -euo pipefail; "
            f"mkdir -p [args[0]]; "
            f"exec > [args[4]] 2>&1; "
            f"jgi_summarize_bam_contig_depths --outputDepth [targets[0]] [depends[0]] && "
            f"metabat2 -i [depends[1]] -a [targets[0]] -o [args[0]]/bin "
            f"-t [args[1]] --minContig [args[2]] [args[3]] && "
            f"touch [targets[1]]"
            f"'",
            depends=[sorted_bam, contigs],
            targets=[depth_file, metabat_done],
            args=[sample_bin_dir, threads, min_contig, metabat_extra_args, log_metabat],
            cores=scheduler_cores,
            mem=mem_mb,
            partition=partition,
            time=1440,
            name=f"binning_{sample}",
        )

        sample_checkm_dir = os.path.join(dir_checkm, sample)
        checkm_stats = os.path.join(sample_checkm_dir, "storage", "bin_stats_ext.tsv")
        checkm_mem = max(mem_mb, 64000)
        log_checkm = os.path.join(sample_checkm_dir, "checkm.log")
        workflow.add_task_gridable(
            f"/bin/bash -c '"
            f"set -euo pipefail; "
            f"mkdir -p [args[2]]; "
            f"checkm lineage_wf -t [args[0]] -x fa [args[1]] [args[2]] > [args[3]] 2>&1"
            f"'",
            depends=metabat_done,
            targets=checkm_stats,
            args=[threads, sample_bin_dir, sample_checkm_dir, log_checkm],
            cores=scheduler_cores,
            mem=checkm_mem,
            partition=partition,
            time=2880,
            name=f"checkm_{sample}",
        )
        all_checkm_stats.append(checkm_stats)

    return all_checkm_stats


def _link_or_copy(source, destination):
    """Create a symlink when possible, otherwise copy the source FASTA."""
    if os.path.exists(destination):
        return
    try:
        os.symlink(os.path.abspath(source), destination)
    except OSError:
        shutil.copy2(source, destination)


def _py_prepare_drep(task):
    """Convert CheckM bin stats into dRep genomeInfo and linked FASTA files."""
    bins_link_dir = task.args[0]
    binning_base_dir = task.args[1]
    checkm_csv_path = task.targets[0].name
    os.makedirs(bins_link_dir, exist_ok=True)
    os.makedirs(os.path.dirname(checkm_csv_path), exist_ok=True)

    all_parsed_data = []
    for checkm_file in task.depends:
        f_path = checkm_file.name
        sample_id = Path(f_path).parents[1].name
        try:
            df = pd.read_csv(f_path, sep="\t", header=None, names=["bin_id", "dict_string"])
            for _, row in df.iterrows():
                stats = ast.literal_eval(row["dict_string"])
                bin_filename = f"{row['bin_id']}.fa"
                original_bin_path = os.path.join(binning_base_dir, sample_id, bin_filename)
                if not os.path.exists(original_bin_path):
                    print(f"Warning: Bin fasta not found, skipping: {original_bin_path}")
                    continue

                new_link_name = f"{sample_id}_{bin_filename}"
                link_path = os.path.join(bins_link_dir, new_link_name)
                _link_or_copy(original_bin_path, link_path)
                all_parsed_data.append(
                    {
                        "genome": new_link_name,
                        "completeness": float(stats["Completeness"]),
                        "contamination": float(stats["Contamination"]),
                    }
                )
        except Exception as exc:
            raise RuntimeError(f"Error parsing CheckM file {f_path}: {exc}") from exc

    if not all_parsed_data:
        raise RuntimeError("No valid bins found for dRep.")

    pd.DataFrame(all_parsed_data).to_csv(checkm_csv_path, index=False)


def run_drep(
    workflow,
    checkm_results,
    drep_work_dir,
    binning_base_dir,
    threads=32,
    scheduler_cores=None,
    mem_mb=200000,
    partition=None,
    completeness=50,
    contamination=5,
    secondary_ani=0.95,
    primary_ani=0.90,
    coverage=0.30,
    extra_args="",
):
    """Create dRep preparation and dereplication tasks."""
    scheduler_cores = scheduler_cores or threads
    bins_link_dir = os.path.join(drep_work_dir, "all_bins_linked")
    drep_info_csv = os.path.join(drep_work_dir, "checkm2drep.csv")

    workflow.add_task(
        _py_prepare_drep,
        depends=checkm_results,
        targets=drep_info_csv,
        args=[bins_link_dir, binning_base_dir],
        name="prepare_drep_inputs",
    )

    drep_out_genomes = os.path.join(drep_work_dir, "dereplicated_genomes")
    drep_done = os.path.join(drep_work_dir, "drep.done")
    log_drep = os.path.join(drep_work_dir, "drep.log")
    workflow.add_task_gridable(
        f"/bin/bash -c '"
        f"set -euo pipefail; "
        f"mkdir -p [args[0]] [args[2]]; "
        f"for item in [args[0]]/*; do "
        f"  [ -e \"$item\" ] || continue; "
        f"  case \"$(basename \"$item\")\" in "
        f"    all_bins_linked|checkm2drep.csv) continue ;; "
        f"    *) rm -rf \"$item\" ;; "
        f"  esac; "
        f"done; "
        f"dRep dereplicate [args[0]] -p [args[1]] -g [args[2]]/*.fa "
        f"--genomeInfo [depends[0]] --completeness [args[3]] --contamination [args[4]] "
        f"-sa [args[5]] -nc [args[6]] -pa [args[7]] [args[8]] > [args[9]] 2>&1; "
        f"touch [targets[0]]"
        f"'",
        depends=drep_info_csv,
        targets=drep_done,
        args=[drep_work_dir, threads, bins_link_dir, completeness, contamination, secondary_ani, coverage, primary_ani, extra_args, log_drep],
        cores=scheduler_cores,
        mem=mem_mb,
        partition=partition,
        time=10080,
        name="drep_run",
    )

    return drep_out_genomes


def run_quantification(
    workflow,
    sgb_dir,
    clean_files_map,
    output_dir,
    threads=16,
    scheduler_cores=None,
    mem_mb=64000,
    partition=None,
    bowtie2_extra_args="--very-sensitive-local",
    coverm_method="relative_abundance",
    coverm_extra_args="",
):
    """Create SGB indexing, read mapping, and CoverM abundance tasks."""
    scheduler_cores = scheduler_cores or threads
    db_dir = os.path.join(output_dir, "database")
    sgb_fasta = os.path.join(db_dir, "SGBs.fna")
    sgb_index_prefix = os.path.join(db_dir, "SGBs_index")
    sgb_index_done = sgb_index_prefix + ".done"
    log_index = os.path.join(db_dir, "bowtie2_build.log")
    drep_done = os.path.join(os.path.dirname(sgb_dir), "drep.done")

    workflow.add_task_gridable(
        f"/bin/bash -c '"
        f"set -euo pipefail; "
        f"mkdir -p [args[0]]; "
        f"cat [args[4]]/*.fa > [targets[0]]; "
        f"bowtie2-build --threads [args[1]] [targets[0]] [args[2]] > [args[3]] 2>&1; "
        f"touch [targets[1]]"
        f"'",
        depends=drep_done,
        targets=[sgb_fasta, sgb_index_done],
        args=[db_dir, threads, sgb_index_prefix, log_index, sgb_dir],
        cores=scheduler_cores,
        mem=mem_mb,
        partition=partition,
        time=1440,
        name="index_sgbs",
    )

    sgb_bams = []
    bam_dir = os.path.join(output_dir, "BAMs")
    for sample, (r1, r2) in clean_files_map.items():
        sgb_bam = os.path.join(bam_dir, f"{sample}.sorted.bam")
        sgb_bam_index = sgb_bam + ".bai"
        log_map = os.path.join(bam_dir, f"{sample}_mapping.log")
        workflow.add_task_gridable(
            f"/bin/bash -c '"
            f"set -euo pipefail; "
            f"mkdir -p [args[0]]; "
            f"exec > [args[3]] 2>&1; "
            f"bowtie2 -p [args[1]] [args[4]] -x [args[2]] -1 [depends[0]] -2 [depends[1]] | "
            f"samtools sort -@ [args[1]] -o [targets[0]] -; "
            f"samtools index [targets[0]] '",
            depends=[r1, r2, sgb_index_done],
            targets=[sgb_bam, sgb_bam_index],
            args=[bam_dir, threads, sgb_index_prefix, log_map, bowtie2_extra_args],
            cores=scheduler_cores,
            mem=mem_mb,
            partition=partition,
            time=2880,
            name=f"map_sgb_{sample}",
        )
        sgb_bams.append(sgb_bam)

    abundance_table = os.path.join(output_dir, "sgb_abundance_matrix.tsv")
    log_coverm = os.path.join(output_dir, "coverm.log")
    drep_done = os.path.join(os.path.dirname(sgb_dir), "drep.done")
    workflow.add_task_gridable(
        f"/bin/bash -c '"
        f"set -euo pipefail; "
        f"mkdir -p [args[0]]; "
        f"coverm genome -m [args[2]] --bam-files [args[0]]/*.sorted.bam "
        f"-f [args[5]]/*.fa -t [args[1]] [args[3]] > [targets[0]] 2> [args[4]]"
        f"'",
        depends=[drep_done] + sgb_bams,
        targets=abundance_table,
        args=[bam_dir, threads, coverm_method, coverm_extra_args, log_coverm, sgb_dir],
        cores=scheduler_cores,
        mem=mem_mb,
        partition=partition,
        time=1440,
        name="coverm_calc",
    )

    return abundance_table


def run_phylophlan_sgb_assignment(
    workflow,
    sgb_dir,
    output_dir,
    database_folder,
    database,
    input_extension=".fa",
    threads=16,
    scheduler_cores=None,
    nproc_io=4,
    mem_mb=64000,
    partition=None,
    clean=False,
    extra_args="",
):
    """Assign dereplicated SGB representatives to PhyloPhlAn SGBs and taxonomy."""
    scheduler_cores = scheduler_cores or threads
    log_file = os.path.join(output_dir, "phylophlan_assign_sgbs.log")
    done_file = os.path.join(output_dir, "phylophlan_assign_sgbs.done")
    clean_flag = "--clean" if clean else ""
    database_arg = f"-d {database}" if database else ""
    drep_done = os.path.join(os.path.dirname(sgb_dir), "drep.done")

    workflow.add_task_gridable(
        f"/bin/bash -c '"
        f"set -euo pipefail; "
        f"mkdir -p [args[0]]; "
        f"phylophlan_assign_sgbs -i [args[1]] -e [args[2]] -o [args[0]] "
        f"--database_folder [args[3]] [args[4]] --nproc_cpu [args[5]] --nproc_io [args[6]] "
        f"[args[7]] [args[8]] > [args[9]] 2>&1; "
        f"touch [targets[0]]"
        f"'",
        depends=drep_done,
        targets=done_file,
        args=[output_dir, sgb_dir, input_extension, database_folder, database_arg, threads, nproc_io, clean_flag, extra_args, log_file],
        cores=scheduler_cores,
        mem=mem_mb,
        partition=partition,
        time=2880,
        name="phylophlan_assign_sgbs",
    )

    return done_file


