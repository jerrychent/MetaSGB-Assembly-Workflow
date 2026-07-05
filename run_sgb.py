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
  --bypass-gtdbtk              Skip GTDB-Tk classification and tree inference.

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


def main():
    """Build and run the SGB assembly workflow."""
    workflow = Workflow(
        version="3.3",
        description="Modular SGB assembly pipeline with configurable resources and bypass support",
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
    workflow.add_argument("bypass-gtdbtk", action="store_true", desc="Skip GTDB-Tk classification and tree inference")

    args = workflow.parse_args()
    suffix_r1, suffix_r2 = _paired_suffixes(args.extension_paired)
    res = ResourceConfig(args.resource_cfg)

    # Define the stable output layout used by all stages.
    dir_clean = os.path.join(args.output, "02_Cleandata")
    dir_assembly = os.path.join(args.output, "03_Assembly")
    dir_binning_base = os.path.join(args.output, "04_Binning")
    dir_checkm_base = os.path.join(args.output, "05_Checkm")
    dir_drep = os.path.join(args.output, "06_dRep_SGBs")
    dir_gtdbtk = os.path.join(args.output, "07_GTDBTk")
    dir_abundance = os.path.join(args.output, "08_Abundance")

    # Create top-level output directories before grid jobs start.
    for directory in [dir_clean, dir_assembly, dir_binning_base, dir_checkm_base, dir_drep, dir_gtdbtk, dir_abundance]:
        os.makedirs(directory, exist_ok=True)

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
            threads=kd_res["cores"],
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
            threads=asm_res["cores"],
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
            threads=bin_res["cores"],
            mem_mb=bin_res["mem"],
            partition=bin_res["partition"],
            min_contig=res.getint("binning", "min_contig", fallback=1500),
            bowtie2_extra_args=res.get("binning", "bowtie2_extra_args", fallback=""),
            metabat_extra_args=res.get("binning", "metabat_extra_args", fallback=""),
        )

    # Step 4: dereplicate all bins into SGB representatives.
    if args.bypass_drep:
        print("Info: Bypassing dRep. Using existing SGB directory...")
        sgb_dir = os.path.join(dir_drep, "dereplicated_genomes")
        if not os.path.exists(sgb_dir) or not glob.glob(os.path.join(sgb_dir, "*.fa")):
            raise FileNotFoundError(f"Bypassed dRep but no SGB FASTA files found in {sgb_dir}")
    else:
        drep_res = res.get_params("drep")
        sgb_dir = tasks.run_drep(
            workflow,
            checkm_results,
            dir_drep,
            binning_base_dir=dir_binning_base,
            threads=drep_res["cores"],
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
            threads=quant_res["cores"],
            mem_mb=quant_res["mem"],
            partition=quant_res["partition"],
            coverm_method=res.get("quantification", "coverm_method", fallback="relative_abundance"),
            coverm_extra_args=res.get("quantification", "coverm_extra_args", fallback=""),
        )

    # Step 6: classify SGB representatives with GTDB-Tk and optionally infer a tree.
    if not args.bypass_gtdbtk:
        gtdb_res = res.get_params("gtdbtk")
        tasks.run_gtdbtk(
            workflow,
            sgb_dir,
            dir_gtdbtk,
            threads=gtdb_res["cores"],
            mem_mb=gtdb_res["mem"],
            partition=gtdb_res["partition"],
            marker_set=res.get("gtdbtk", "marker_set", fallback="bac120"),
            skip_ani_screen=res.getboolean("gtdbtk", "skip_ani_screen", fallback=True),
            infer_tree=res.getboolean("gtdbtk", "infer_tree", fallback=True),
            infer_mem_mb=res.getint("gtdbtk", "infer_memory_mb", fallback=128000),
            extra_args=res.get("gtdbtk", "extra_args", fallback=""),
        )

    workflow.go()


if __name__ == "__main__":
    main()
