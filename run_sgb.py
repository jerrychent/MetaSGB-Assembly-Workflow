import glob
import os
import sys

HELP_TEXT = """MetaSGB Assembly Workflow

Usage:
  python run_sgb.py --input RAW_READ_DIR --output OUTPUT_DIR [options]

Required anadama2 arguments:
  --input PATH                 Directory containing raw paired-end reads.
  --output PATH                Output directory for all workflow results.

Workflow options:
  --resource-cfg PATH          Resource and tool config file.
                               Default: config/resources.cfg
  --extension-paired R1,R2     Paired read suffixes separated by a comma.
                               Default: _1.fastq.gz,_2.fastq.gz
  --kneaddata-db PATH          KneadData database path for host/contaminant removal.
                               Default: /public/databases/kneaddata_db

Bypass options:
  --bypass-kneaddata           Use existing clean reads in 02_Cleandata.
  --bypass-assembly            Use existing contigs in 03_Assembly.
  --bypass-binning             Use existing CheckM outputs in 05_Checkm.
  --bypass-drep                Use existing SGB FASTA files in 06_dRep_SGBs/dereplicated_genomes.
  --bypass-quant               Skip CoverM abundance quantification.
  --bypass-phylophlan-sgb      Skip PhyloPhlAn SGB and taxonomy assignment.

Environment assumption:
  Activate the integrated workflow environment before running this script. Tool-specific
  env prefix arguments are intentionally not exposed.

Examples:
  python run_sgb.py --input raw_reads --output sgb_out --kneaddata-db /db/human
  python run_sgb.py --input raw_reads --output sgb_out --bypass-kneaddata --bypass-assembly
"""


# Print help before importing anadama2 so help works in lightweight environments.
if "--help" in sys.argv[1:] or "-h" in sys.argv[1:]:
    print(HELP_TEXT)
    sys.exit(0)

from anadama2 import Workflow

from lib.config_loader import ResourceConfig
from lib.slurm_patch import patch_slurm_cpus_per_task
import lib.sgb_tasks as tasks


def _paired_suffixes(value):
    """Parse the R1/R2 suffix pair from the CLI string."""
    suffixes = [item.strip() for item in value.split(",") if item.strip()]
    if len(suffixes) != 2:
        raise ValueError(
            "Error: --extension-paired must contain two suffixes separated by a comma "
            "(e.g. '_1.fastq.gz,_2.fastq.gz')."
        )
    return suffixes


def _find_clean_pairs(clean_dir):
    """Find existing KneadData paired reads for bypass mode."""
    clean_files_info = []
    patterns = ["*_paired_1.fastq", "*_paired_1.fastq.gz"]
    seen = set()

    for pattern in patterns:
        for r1 in sorted(glob.glob(os.path.join(clean_dir, "**", pattern), recursive=True)):
            if r1 in seen:
                continue
            seen.add(r1)

            if r1.endswith("_paired_1.fastq.gz"):
                sample = os.path.basename(r1).replace("_paired_1.fastq.gz", "")
                r2 = r1.replace("_paired_1.fastq.gz", "_paired_2.fastq.gz")
            else:
                sample = os.path.basename(r1).replace("_paired_1.fastq", "")
                r2 = r1.replace("_paired_1.fastq", "_paired_2.fastq")

            if os.path.exists(r2):
                clean_files_info.append((sample, r1, r2))
            else:
                print(f"Warning: Missing paired clean R2 for {sample}: {r2}")

    return clean_files_info


def _find_assemblies(assembly_dir, clean_files_info):
    """Find existing contig files for assembly bypass mode."""
    assembly_info = []
    for sample, _, _ in clean_files_info:
        contigs = os.path.join(assembly_dir, sample, "contigs.fasta")
        if os.path.exists(contigs):
            assembly_info.append((sample, contigs))
        else:
            print(f"Warning: Contigs for {sample} not found: {contigs}")
    return assembly_info


def _remove_file_if_exists(path):
    """Remove a file if it exists."""
    try:
        os.remove(path)
    except FileNotFoundError:
        pass


def _touch_file(path):
    """Create an empty file or update its mtime."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "a", encoding="utf-8"):
        os.utime(path, None)


def _latest_nonlog_mtime(directory):
    """Return the newest mtime from result files under a directory."""
    if not os.path.isdir(directory):
        return None
    latest = None
    for root, _, files in os.walk(directory):
        for name in files:
            if name.endswith(".log") or name.endswith(".done"):
                continue
            path = os.path.join(root, name)
            try:
                mtime = os.path.getmtime(path)
            except FileNotFoundError:
                continue
            latest = mtime if latest is None else max(latest, mtime)
    return latest


def _mark_done_from_existing_outputs(done_file, output_dir):
    """Create a done marker using the newest existing result timestamp."""
    result_mtime = _latest_nonlog_mtime(output_dir)
    os.makedirs(os.path.dirname(done_file), exist_ok=True)
    with open(done_file, "a", encoding="utf-8"):
        pass
    if result_mtime is not None:
        os.utime(done_file, (result_mtime, result_mtime))


def _directory_has_fastas(directory):
    """Return True when a directory contains at least one FASTA-like file."""
    if not os.path.isdir(directory):
        return False
    patterns = ("*.fa", "*.fna", "*.fasta", "*.fa.gz", "*.fna.gz", "*.fasta.gz")
    return any(glob.glob(os.path.join(directory, pattern)) for pattern in patterns)


def _directory_has_nonlog_entries(directory):
    """Return True when a directory contains files other than logs and done markers."""
    if not os.path.isdir(directory):
        return False
    for name in os.listdir(directory):
        if name.endswith(".log") or name.endswith(".done"):
            continue
        return True
    return False


def _drep_is_complete(drep_dir):
    """Detect a completed dRep run from legacy or current outputs."""
    drep_done = os.path.join(drep_dir, "drep.done")
    drep_csv = os.path.join(drep_dir, "checkm2drep.csv")
    sgb_dir = os.path.join(drep_dir, "dereplicated_genomes")
    return _directory_has_fastas(sgb_dir) and (os.path.exists(drep_done) or os.path.exists(drep_csv))


def _phylophlan_is_complete(output_dir):
    """Detect a completed PhyloPhlAn assignment run."""
    done_file = os.path.join(output_dir, "phylophlan_assign_sgbs.done")
    return os.path.exists(done_file) or _directory_has_nonlog_entries(output_dir)


def main():
    """Build and run the SGB assembly workflow."""
    patch_slurm_cpus_per_task()
    workflow = Workflow(
        version="3.4",
        description="Modular SGB assembly pipeline with PhyloPhlAn SGB assignment",
    )

    # Register workflow-level arguments.
    workflow.add_argument("resource-cfg", desc="Path to resource config file", default="config/resources.cfg")
    workflow.add_argument(
        "extension-paired",
        desc="Suffixes for paired reads, comma separated (e.g. '_1.fastq.gz,_2.fastq.gz')",
        default="_1.fastq.gz,_2.fastq.gz",
    )
    workflow.add_argument(
        "kneaddata-db",
        desc="Path to the KneadData database for the target host/contaminant species",
        default="/public/databases/kneaddata_db",
    )

    # Bypass switches allow resuming from completed stages.
    workflow.add_argument("bypass-kneaddata", action="store_true", desc="Use existing clean reads in 02_Cleandata")
    workflow.add_argument("bypass-assembly", action="store_true", desc="Use existing contigs in 03_Assembly")
    workflow.add_argument("bypass-binning", action="store_true", desc="Use existing CheckM outputs in 05_Checkm")
    workflow.add_argument("bypass-drep", action="store_true", desc="Use existing SGB FASTA files in 06_dRep_SGBs/dereplicated_genomes")
    workflow.add_argument("bypass-quant", action="store_true", desc="Skip abundance quantification")
    workflow.add_argument("bypass-phylophlan-sgb", action="store_true", desc="Skip PhyloPhlAn SGB and taxonomy assignment")

    args = workflow.parse_args()
    suffix_r1, suffix_r2 = _paired_suffixes(args.extension_paired)
    res = ResourceConfig(args.resource_cfg)

    # Start each workflow run with a fresh anadama2 log.
    _remove_file_if_exists(os.path.join(args.output, "anadama.log"))
    _remove_file_if_exists(os.path.join(os.getcwd(), "anadama.log"))

    # Define the stable output layout used by all stages.
    dir_clean = os.path.join(args.output, "02_Cleandata")
    dir_assembly = os.path.join(args.output, "03_Assembly")
    dir_binning_base = os.path.join(args.output, "04_Binning")
    dir_checkm_base = os.path.join(args.output, "05_Checkm")
    dir_drep = os.path.join(args.output, "06_dRep_SGBs")
    dir_abundance = os.path.join(args.output, "08_Abundance")
    dir_phylophlan = os.path.join(args.output, "09_PhyloPhlAn_SGB_Assignment")

    # Create top-level output directories before grid jobs start.
    for directory in [dir_clean, dir_assembly, dir_binning_base, dir_checkm_base, dir_drep, dir_abundance, dir_phylophlan]:
        os.makedirs(directory, exist_ok=True)

    # Migrate legacy completed results to stable done markers so anadama2 can skip them.
    if _drep_is_complete(dir_drep):
        _mark_done_from_existing_outputs(os.path.join(dir_drep, "drep.done"), os.path.join(dir_drep, "dereplicated_genomes"))
    if _phylophlan_is_complete(dir_phylophlan):
        _mark_done_from_existing_outputs(os.path.join(dir_phylophlan, "phylophlan_assign_sgbs.done"), dir_phylophlan)

    # Step 1: clean and decontaminate paired-end reads with KneadData.
    if args.bypass_kneaddata:
        print("Info: Bypassing KneadData. Searching for existing clean files...")
        clean_files_info = _find_clean_pairs(dir_clean)
        if not clean_files_info:
            raise FileNotFoundError(f"Bypassed KneadData but no paired clean reads found under {dir_clean}")
    else:
        input_files = workflow.get_input_files(extension=suffix_r1)
        if not input_files:
            raise FileNotFoundError(f"No input files found with R1 suffix: {suffix_r1}")

        kd_res = res.get_params("kneaddata")
        clean_files_info = tasks.run_kneaddata(
            workflow,
            input_files,
            dir_clean,
            db_path=args.kneaddata_db,
            threads=kd_res["threads"],
            scheduler_cores=kd_res["cores"],
            mem_mb=kd_res["mem"],
            partition=kd_res["partition"],
            suffix_r1=suffix_r1,
            suffix_r2=suffix_r2,
            sequencer_source=res.get("kneaddata", "sequencer_source", fallback="TruSeq3"),
            run_fastqc_start=res.getboolean("kneaddata", "run_fastqc_start", fallback=True),
            run_fastqc_end=res.getboolean("kneaddata", "run_fastqc_end", fallback=True),
            extra_args=res.get("kneaddata", "extra_args", fallback=""),
        )

    # Step 2: assemble clean reads into per-sample contigs.
    if args.bypass_assembly:
        print("Info: Bypassing Assembly. Searching for existing contigs...")
        assembly_info = _find_assemblies(dir_assembly, clean_files_info)
        if not assembly_info:
            raise FileNotFoundError(f"Bypassed Assembly but no contigs found under {dir_assembly}")
    else:
        asm_res = res.get_params("assembly")
        assembly_info = tasks.run_assembly(
            workflow,
            clean_files_info,
            dir_assembly,
            threads=asm_res["threads"],
            scheduler_cores=asm_res["cores"],
            mem_mb=asm_res["mem"],
            partition=asm_res["partition"],
            extra_args=res.get("assembly", "extra_args", fallback=""),
        )

    # Step 3: map reads to contigs, bin assemblies, and run CheckM.
    clean_files_map = {sample: (r1, r2) for sample, r1, r2 in clean_files_info}
    if args.bypass_binning:
        print("Info: Bypassing Binning. Searching for existing CheckM stats...")
        checkm_results = []
        for sample, _ in assembly_info:
            stats_file = os.path.join(dir_checkm_base, sample, "storage", "bin_stats_ext.tsv")
            if os.path.exists(stats_file):
                checkm_results.append(stats_file)
            else:
                print(f"Warning: CheckM stats for {sample} not found: {stats_file}")
        if not checkm_results:
            raise FileNotFoundError(f"Bypassed Binning but no CheckM results found under {dir_checkm_base}")
    else:
        bin_res = res.get_params("binning")
        checkm_results = tasks.run_binning_workflow(
            workflow,
            assembly_info,
            clean_files_map,
            args.output,
            threads=bin_res["threads"],
            scheduler_cores=bin_res["cores"],
            mem_mb=bin_res["mem"],
            partition=bin_res["partition"],
            min_contig=res.getint("binning", "min_contig", fallback=1500),
            bowtie2_extra_args=res.get("binning", "bowtie2_extra_args", fallback="--very-sensitive-local"),
            metabat_extra_args=res.get("binning", "metabat_extra_args", fallback=""),
        )

    # Step 4: dereplicate all bins into SGB representatives.
    if args.bypass_drep:
        if _drep_is_complete(dir_drep):
            print("Info: Bypassing dRep. Using existing dRep results.")
            sgb_dir = os.path.join(dir_drep, "dereplicated_genomes")
            _mark_done_from_existing_outputs(os.path.join(dir_drep, "drep.done"), sgb_dir)
        else:
            raise FileNotFoundError(f"Bypassed dRep but no complete results found under {dir_drep}")
    else:
        drep_res = res.get_params("drep")
        sgb_dir = tasks.run_drep(
            workflow,
            checkm_results,
            dir_drep,
            binning_base_dir=dir_binning_base,
            threads=drep_res["threads"],
            scheduler_cores=drep_res["cores"],
            mem_mb=drep_res["mem"],
            partition=drep_res["partition"],
            completeness=res.getfloat("drep", "completeness", fallback=50),
            contamination=res.getfloat("drep", "contamination", fallback=5),
            secondary_ani=res.getfloat("drep", "secondary_ani", fallback=0.95),
            primary_ani=res.getfloat("drep", "primary_ani", fallback=0.90),
            coverage=res.getfloat("drep", "coverage", fallback=0.30),
            extra_args=res.get("drep", "extra_args", fallback=""),
        )

    # Step 5: quantify SGB abundance across all samples.
    if not args.bypass_quant:
        quant_res = res.get_params("quantification", fallback_step="binning")
        tasks.run_quantification(
            workflow,
            sgb_dir,
            clean_files_map,
            dir_abundance,
            threads=quant_res["threads"],
            scheduler_cores=quant_res["cores"],
            mem_mb=quant_res["mem"],
            partition=quant_res["partition"],
            bowtie2_extra_args=res.get("quantification", "bowtie2_extra_args", fallback="--very-sensitive-local"),
            coverm_method=res.get("quantification", "coverm_method", fallback="relative_abundance"),
            coverm_extra_args=res.get("quantification", "coverm_extra_args", fallback=""),
        )

    # Step 6: assign SGB and taxonomy with PhyloPhlAn.
    if args.bypass_phylophlan_sgb:
        if _phylophlan_is_complete(dir_phylophlan):
            print("Info: Bypassing PhyloPhlAn. Using existing SGB assignment results.")
            _mark_done_from_existing_outputs(os.path.join(dir_phylophlan, "phylophlan_assign_sgbs.done"), dir_phylophlan)
        else:
            raise FileNotFoundError(f"Bypassed PhyloPhlAn but no complete results found under {dir_phylophlan}")
    else:
        phy_res = res.get_params("phylophlan_sgb")
        tasks.run_phylophlan_sgb_assignment(
            workflow,
            sgb_dir,
            dir_phylophlan,
            database_folder=res.get("phylophlan_sgb", "database_folder", fallback="phylophlan_databases"),
            database=res.get("phylophlan_sgb", "database", fallback=""),
            input_extension=res.get("phylophlan_sgb", "input_extension", fallback=".fa"),
            threads=phy_res["threads"],
            scheduler_cores=phy_res["cores"],
            nproc_io=res.getint("phylophlan_sgb", "nproc_io", fallback=4),
            mem_mb=phy_res["mem"],
            partition=phy_res["partition"],
            clean=res.getboolean("phylophlan_sgb", "clean", fallback=False),
            extra_args=res.get("phylophlan_sgb", "extra_args", fallback=""),
        )

    workflow.go()


if __name__ == "__main__":
    main()
