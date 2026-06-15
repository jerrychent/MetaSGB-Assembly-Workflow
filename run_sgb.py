import os
import glob
from anadama2 import Workflow
from lib.config_loader import ResourceConfig
import lib.sgb_tasks as tasks

def main():
    # --- 1. 初始化工作流 ---
    workflow = Workflow(
        version="3.0", 
        description="Modular SGB Pipeline with Bypass Support"
    )

    # --- 2. 定义命令行参数 ---
    # 2.1 核心路径参数
    workflow.add_argument("resource-cfg", desc="Path to resource config file", default="config/resources.cfg")

    workflow.add_argument("extension-paired", 
                          desc="Suffixes for paired reads, comma separated (e.g. '_1.fastq.gz,_2.fastq.gz')", 
                          default="_1.fastq.gz,_2.fastq.gz")
    
    # 2.2 数据库与软件路径
    workflow.add_argument("kneaddata-db", desc="Path to KneadData Database", 
                          default="/public/databases/kneaddata_db")
#    workflow.add_argument("fastqc-path", desc="Path to FastQC executable", 
#                          default="fastqc") 
    #workflow.add_argument("trimmomatic-path", desc="Path to Trimmomatic JAR directory", required=True)
    
    # 2.3 环境设置
    workflow.add_argument("env-kneaddata", default="", desc="Env string for KneadData")
    workflow.add_argument("env-spades",    default="", desc="Env string for MetaSPAdes")
    workflow.add_argument("env-bowtie",    default="", desc="Env string for Bowtie2/Samtools")
    workflow.add_argument("env-binning",   default="", desc="Env string for MetaBAT2/CheckM")
    workflow.add_argument("env-drep",      default="", desc="Env string for dRep")
    workflow.add_argument("env-coverm",    default="", desc="Env string for CoverM")
    workflow.add_argument("env-gtdbtk",    default="", desc="Env string for GTDB-Tk")

    # 2.4 Bypass 参数
    workflow.add_argument("bypass-kneaddata", action="store_true", desc="Skip KneadData step")
    workflow.add_argument("bypass-assembly",  action="store_true", desc="Skip Assembly step")
    workflow.add_argument("bypass-binning",   action="store_true", desc="Skip Binning & CheckM step")
    workflow.add_argument("bypass-drep",      action="store_true", desc="Skip dRep step")
    workflow.add_argument("bypass-quant",     action="store_true", desc="Skip Quantification step")

    args = workflow.parse_args()

    try:
        suffix_r1, suffix_r2 = [s.strip() for s in args.extension_paired.split(',')]
    except ValueError:
        raise ValueError("Error: --extension-paired must contain two suffixes separated by a comma (e.g. '_1.fastq.gz,_2.fastq.gz')")
    
    # --- 3. 加载配置文件 ---

    res = ResourceConfig(args.resource_cfg)

    # --- 4. 准备路径 ---
    # 定义各步骤输出目录
    dir_clean = os.path.join(args.output, "02_Cleandata")
    dir_assembly = os.path.join(args.output, "03_Assembly")
    dir_binning_base = os.path.join(args.output, "04_Binning") 
    dir_checkm_base = os.path.join(args.output, "05_Checkm")
    dir_drep = os.path.join(args.output, "06_dRep_SGBs")
    dir_gtdbtk = os.path.join(args.output, "07_GTDBTk")
    dir_abundance = os.path.join(args.output, "08_Abundance")

    dirs_to_create = [
        dir_clean, 
        dir_assembly, 
        dir_binning_base, 
        dir_checkm_base, 
        dir_drep, 
        dir_gtdbtk, 
        dir_abundance
    ]
    
    print("-" * 50)
    print("Info: Pre-creating output directories to prevent race conditions...")
    for d in dirs_to_create:
        if not os.path.exists(d):
            try:
                os.makedirs(d, exist_ok=True)
                print(f"Created: {d}")
            except Exception as e:
                print(f"Error creating directory {d}: {e}")
                # 如果创建失败，通常是权限问题，提前报错退出比跑到一半报错好
                raise e
    print("-" * 50)

    # ==========================================================================
    # 流程执行阶段 
    # ==========================================================================

    # --- Step 1: KneadData 质控 ---
    clean_files_info = []
    
    if not args.bypass_kneaddata:
        input_files = workflow.get_input_files(extension=suffix_r1)
        kd_res = res.get_params("kneaddata")
        
        clean_files_info = tasks.run_kneaddata(
            workflow, 
            input_files, 
            dir_clean,
            db_path=args.kneaddata_db,
            threads=kd_res['cores'],
            mem_mb=kd_res['mem'],
            partition=kd_res['partition'],
            suffix_r1=suffix_r1, 
            suffix_r2=suffix_r2,
            env=args.env_kneaddata
        )
    else:
        print("Info: Bypassing KneadData. Searching for existing clean files...")
        found_r1 = glob.glob(os.path.join(dir_clean, "*_paired_1.fastq"))
        for r1 in found_r1:
            sample = os.path.basename(r1).replace("_paired_1.fastq", "")
            r2 = r1.replace("_paired_1.fastq", "_paired_2.fastq")
            if os.path.exists(r2):
                clean_files_info.append((sample, r1, r2))
        
        if not clean_files_info:
            raise FileNotFoundError(f"Bypassed KneadData but no files found in {dir_clean}")

    # --- Step 2: MetaSPAdes 组装 ---
    assembly_info = []
    
    if not args.bypass_assembly:
        asm_res = res.get_params("assembly")
        assembly_info = tasks.run_assembly(
            workflow,
            clean_files_info,
            dir_assembly,
            threads=asm_res['cores'],
            mem_mb=asm_res['mem'],
            partition=asm_res['partition'],
            env=args.env_spades
        )
    else:
        print("Info: Bypassing Assembly. Searching for existing contigs...")
        for sample, _, _ in clean_files_info:
            contigs = os.path.join(dir_assembly, sample, "contigs.fasta")
            if os.path.exists(contigs):
                assembly_info.append((sample, contigs))
            else:
                print(f"Warning: Contigs for {sample} not found in {dir_assembly}")

        if not assembly_info:
            raise FileNotFoundError(f"Bypassed Assembly but no contigs found in {dir_assembly}")

    # --- Step 3 & 4: Binning ---
    checkm_results = []
    clean_files_map = {item[0]: (item[1], item[2]) for item in clean_files_info}
    
    if not args.bypass_binning:
        bin_res = res.get_params("binning")
        checkm_results = tasks.run_binning_workflow(
            workflow,
            assembly_info,
            clean_files_map,
            args.output,
            threads=bin_res['cores'],
            mem_mb=bin_res['mem'],
            partition=bin_res['partition'],
            env_bowtie=args.env_bowtie,
            env_binning=args.env_binning
        )
    else:
        print("Info: Bypassing Binning. Searching for existing CheckM stats...")
        for sample, _ in assembly_info:
            stats_file = os.path.join(dir_checkm_base, sample, "storage", "bin_stats_ext.tsv")
            if os.path.exists(stats_file):
                checkm_results.append(stats_file) 
        if not checkm_results:
             raise FileNotFoundError(f"Bypassed Binning but no CheckM results found in {dir_checkm_base}")

    # --- Step 5: dRep ---
    sgb_dir = ""
    if not args.bypass_drep:
        drep_res = res.get_params("drep")
        binning_base_dir = dir_binning_base
        sgb_dir = tasks.run_drep(
            workflow,
            checkm_results, 
            dir_drep,
            binning_base_dir=binning_base_dir,
            threads=drep_res['cores'],
            mem_mb=drep_res['mem'],
            partition=drep_res['partition'],
            env=args.env_drep
        )
    else:
        print("Info: Bypassing dRep. Using existing SGB directory...")
        sgb_dir = os.path.join(dir_drep, "dereplicated_genomes")
        if not os.path.exists(sgb_dir) or not glob.glob(os.path.join(sgb_dir, "*.fa")):
             raise FileNotFoundError(f"Bypassed dRep but no SGBs found in {sgb_dir}")

    # --- Step 6: Quantification ---
    if not args.bypass_quant:
        quant_res = res.get_params("binning") 
        tasks.run_quantification(
            workflow,
            sgb_dir,
            clean_files_map,
            dir_abundance,
            threads=quant_res['cores'],
            mem_mb=64000,
            partition=quant_res['partition'],
            env_bowtie=args.env_bowtie,
            env_coverm=args.env_coverm
        )

    # --- Step 7: GTDB-Tk ---
    gtdb_res = res.get_params("gtdbtk")
    tasks.run_gtdbtk(
        workflow,
        sgb_dir,
        dir_gtdbtk,
        threads=gtdb_res['cores'],
        mem_mb=gtdb_res['mem'],
        partition=gtdb_res['partition'],
        env=args.env_gtdbtk
    )

    workflow.go()

if __name__ == "__main__":
    main()