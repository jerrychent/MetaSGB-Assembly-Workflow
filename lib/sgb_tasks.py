import ast
import os
import shutil
from pathlib import Path

import pandas as pd
from anadama2.tracked import TrackedDirectory


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
    mem_mb=10240,
    partition=None,
    sequencer_source="TruSeq3",
    run_fastqc_start=True,
    run_fastqc_end=True,
    extra_args="",
):
    """Create one KneadData task per paired-end sample."""
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
            f"kneaddata --input1 [depends[0]] --input2 [depends[1]] "
            f"--output [args[0]] --output-prefix [args[1]] "
            f"-t [args[2]] -db [args[3]] --bowtie2-options=\"-p [args[2]]\" "
            f"--sequencer-source [args[4]] [args[5]] [args[6]] "
            f"> [args[7]] 2>&1"
            f"'"
        )

        workflow.add_task_gridable(
            cmd,
            depends=[r1, r2],
            targets=[clean_r1, clean_r2],
            args=[sample_out_dir, sample, threads, db_path, sequencer_source, fastqc_flags, extra_args, log_file],
            cores=threads,
            mem=mem_mb,
            partition=partition,
            time=2880,
            name=f"knead_{sample}",
        )
        clean_files_info.append((sample, clean_r1, clean_r2))

    return clean_files_info


def run_assembly(workflow, clean_files_info, output_base_dir, threads=16, mem_mb=200000, partition=None, extra_args=""):
    """Create one metaSPAdes assembly task per sample."""
    assembly_info = []

    for sample, r1, r2 in clean_files_info:
        sample_out_dir = os.path.join(output_base_dir, sample)
        contigs = os.path.join(sample_out_dir, "contigs.fasta")
        log_file = os.path.join(sample_out_dir, "spades_console.log")

        cmd = (
            f"/bin/bash -c '"
            f"set -euo pipefail; "
            f"mkdir -p [args[0]]; "
            f"metaspades.py -1 [depends[0]] -2 [depends[1]] "
            f"-o [args[0]] -t [args[1]] --memory [args[2]] [args[3]] "
            f"> [args[4]] 2>&1"
            f"'"
        )

        workflow.add_task_gridable(
            cmd,
            depends=[r1, r2],
            targets=contigs,
            args=[sample_out_dir, threads, int(mem_mb / 1024), extra_args, log_file],
            cores=threads,
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
    mem_mb=32000,
    partition=None,
    min_contig=1500,
    bowtie2_extra_args="",
    metabat_extra_args="",
):
    """Create Bowtie2, MetaBAT2, and CheckM tasks for each assembly."""
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
            cores=threads,
            mem=mem_mb,
            partition=partition,
            time=1440,
            name=f"index_{sample}",
        )

        # Pipefail makes Bowtie2 failures propagate through samtools sort.
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
            cores=threads,
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
            f"jgi_summarize_bam_contig_depths --outputDepth [targets[0]] [depends[0]] && "
            f"metabat2 -i [depends[1]] -a [targets[0]] -o [args[0]]/bin "
            f"-t [args[1]] --minContig [args[2]] [args[3]] && "
            f"touch [targets[1]]"
            f"' > [args[4]] 2>&1",
            depends=[sorted_bam, contigs],
            targets=[depth_file, metabat_done],
            args=[sample_bin_dir, threads, min_contig, metabat_extra_args, log_metabat],
            cores=threads,
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
            cores=threads,
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
    log_drep = os.path.join(drep_work_dir, "drep.log")
    workflow.add_task_gridable(
        f"/bin/bash -c '"
        f"set -euo pipefail; "
        f"mkdir -p [args[0]]; "
        f"dRep dereplicate [args[0]] -p [args[1]] -g [args[2]]/*.fa "
        f"--genomeInfo [depends[0]] --completeness [args[3]] --contamination [args[4]] "
        f"-sa [args[5]] -nc [args[6]] -pa [args[7]] [args[8]] > [args[9]] 2>&1"
        f"'",
        depends=drep_info_csv,
        targets=TrackedDirectory(drep_out_genomes),
        args=[drep_work_dir, threads, bins_link_dir, completeness, contamination, secondary_ani, coverage, primary_ani, extra_args, log_drep],
        cores=threads,
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
    mem_mb=64000,
    partition=None,
    coverm_method="relative_abundance",
    coverm_extra_args="",
):
    """Create SGB indexing, read mapping, and CoverM abundance tasks."""
    db_dir = os.path.join(output_dir, "database")
    sgb_fasta = os.path.join(db_dir, "SGBs.fna")
    sgb_index_prefix = os.path.join(db_dir, "SGBs_index")
    sgb_index_done = sgb_index_prefix + ".done"
    log_index = os.path.join(db_dir, "bowtie2_build.log")

    workflow.add_task_gridable(
        f"/bin/bash -c '"
        f"set -euo pipefail; "
        f"mkdir -p [args[0]]; "
        f"cat [depends[0]]/*.fa > [targets[0]]; "
        f"bowtie2-build --threads [args[1]] [targets[0]] [args[2]] > [args[3]] 2>&1; "
        f"touch [targets[1]]"
        f"'",
        depends=TrackedDirectory(sgb_dir),
        targets=[sgb_fasta, sgb_index_done],
        args=[db_dir, threads, sgb_index_prefix, log_index],
        cores=threads,
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
            f"bowtie2 -p [args[1]] -x [args[2]] -1 [depends[0]] -2 [depends[1]] | "
            f"samtools sort -@ [args[1]] -o [targets[0]] -; "
            f"samtools index [targets[0]]"
            f"' > [args[3]] 2>&1",
            depends=[r1, r2, sgb_index_done],
            targets=[sgb_bam, sgb_bam_index],
            args=[bam_dir, threads, sgb_index_prefix, log_map],
            cores=threads,
            mem=mem_mb,
            partition=partition,
            time=2880,
            name=f"map_sgb_{sample}",
        )
        sgb_bams.append(sgb_bam)

    abundance_table = os.path.join(output_dir, "sgb_abundance_matrix.tsv")
    log_coverm = os.path.join(output_dir, "coverm.log")
    workflow.add_task_gridable(
        f"/bin/bash -c '"
        f"set -euo pipefail; "
        f"mkdir -p [args[0]]; "
        f"coverm genome -m [args[2]] --bam-files [args[0]]/*.sorted.bam "
        f"-f [depends[0]]/*.fa -t [args[1]] [args[3]] > [targets[0]] 2> [args[4]]"
        f"'",
        depends=[TrackedDirectory(sgb_dir)] + sgb_bams,
        targets=abundance_table,
        args=[bam_dir, threads, coverm_method, coverm_extra_args, log_coverm],
        cores=threads,
        mem=mem_mb,
        partition=partition,
        time=1440,
        name="coverm_calc",
    )

    return abundance_table


def run_gtdbtk(
    workflow,
    sgb_dir,
    output_dir,
    threads=32,
    mem_mb=500000,
    partition=None,
    marker_set="bac120",
    skip_ani_screen=True,
    infer_tree=True,
    infer_mem_mb=128000,
    extra_args="",
):
    """Create GTDB-Tk classification and optional tree inference tasks."""
    gtdb_summary = os.path.join(output_dir, f"gtdbtk.{marker_set}.summary.tsv")
    msa_file = os.path.join(output_dir, "align", f"gtdbtk.{marker_set}.user_msa.fasta")
    log_classify = os.path.join(output_dir, "gtdbtk_classify.log")
    skip_ani_flag = "--skip_ani_screen" if skip_ani_screen else ""

    workflow.add_task_gridable(
        f"/bin/bash -c '"
        f"set -euo pipefail; "
        f"mkdir -p [args[0]]; "
        f"gtdbtk classify_wf --genome_dir [depends[0]] --out_dir [args[0]] "
        f"--extension fa --cpus [args[1]] [args[2]] [args[3]] > [args[4]] 2>&1"
        f"'",
        depends=TrackedDirectory(sgb_dir),
        targets=[gtdb_summary, msa_file],
        args=[output_dir, threads, skip_ani_flag, extra_args, log_classify],
        cores=threads,
        mem=mem_mb,
        partition=partition,
        time=10080,
        name="gtdbtk_classify",
    )

    if not infer_tree:
        return gtdb_summary

    tree_out_dir = os.path.join(output_dir, "classify")
    log_infer = os.path.join(output_dir, "gtdbtk_infer.log")
    workflow.add_task_gridable(
        f"/bin/bash -c '"
        f"set -euo pipefail; "
        f"mkdir -p [args[0]]; "
        f"gtdbtk infer --out_dir [args[0]] --cpus [args[1]] --msa_file [depends[0]] > [args[2]] 2>&1"
        f"'",
        depends=msa_file,
        targets=TrackedDirectory(tree_out_dir),
        args=[tree_out_dir, threads, log_infer],
        cores=threads,
        mem=infer_mem_mb,
        partition=partition,
        time=2880,
        name="gtdbtk_infer",
    )

    return gtdb_summary

