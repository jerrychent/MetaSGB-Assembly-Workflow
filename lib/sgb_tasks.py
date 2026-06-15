import os
import shutil
import pandas as pd
import ast
from anadama2.tracked import TrackedDirectory, TrackedExecutable

"""
SGB Pipeline Tasks Library (lib/sgb_tasks.py)
包含 SGB 分析流程中所有具体的任务实现函数。

"""

def run_kneaddata(workflow, input_files, output_dir, db_path, 
                  suffix_r1="_1.fastq.gz", suffix_r2="_2.fastq.gz",
                  threads=4, mem_mb=10240, partition=None, env=""):
    """
    [Step 1] KneadData 质控
    """
    clean_files_info = []
    
    for r1 in input_files:
        dirname = os.path.dirname(r1)
        basename = os.path.basename(r1)
        
        # 动态后缀匹配
        if not basename.endswith(suffix_r1):
            continue
            
        sample = basename.rsplit(suffix_r1, 1)[0]
        basename_r2 = basename.replace(suffix_r1, suffix_r2)
        r2 = os.path.join(dirname, basename_r2)
        
        if not os.path.exists(r2):
            raise FileNotFoundError(f"R2 file not found: {r2}")

        sample_out_dir = os.path.join(output_dir, sample)

        clean_r1 = os.path.join(sample_out_dir, f"{sample}_paired_1.fastq")
        clean_r2 = os.path.join(sample_out_dir, f"{sample}_paired_2.fastq")
        

        log_file = os.path.join(sample_out_dir, f"{sample}_kneaddata.log")

        cmd = (
            f"mkdir -p [args[0]] && " 
            f"{env} kneaddata --input1 [depends[0]] --input2 [depends[1]] "
            f"--output [args[0]] --output-prefix [args[1]] "
            f"-t [args[2]] -db [args[3]] "
            f"--bowtie2-options='-p [args[2]]' --sequencer-source TruSeq3 "
            f"--run-fastqc-start --run-fastqc-end "
            f"> {log_file} 2>&1"
        )

        workflow.add_task_gridable(
            cmd,
            depends=[r1, r2],
            targets=[clean_r1, clean_r2],
            # [关键修改] args[0] 传入 sample_out_dir 而不是大的 output_dir
            args=[sample_out_dir, sample, threads, db_path, mem_mb],
            cores=threads,
            mem=mem_mb,
            partition=partition,
            time=2880,
            name=f"knead_{sample}"
        )
        
        clean_files_info.append((sample, clean_r1, clean_r2))
        
    return clean_files_info

def run_assembly(workflow, clean_files_info, output_base_dir, 
                 threads=16, mem_mb=200000, partition=None, env=""):
    """
    [Step 2] MetaSPAdes 组装
    """
    assembly_info = []
    
    for sample, r1, r2 in clean_files_info:
        sample_out_dir = os.path.join(output_base_dir, sample)
        contigs = os.path.join(sample_out_dir, "contigs.fasta")
        
        log_file = os.path.join(sample_out_dir, "spades_console.log")

        cmd = (
            f"mkdir -p [args[0]] && "
            f"{env} metaspades.py -1 [depends[0]] -2 [depends[1]] "
            f"-o [args[0]] -t [args[1]] --memory [args[2]] "
            f"> {log_file} 2>&1"
        )
        
        workflow.add_task_gridable(
            cmd,
            depends=[r1, r2],
            targets=contigs,
            args=[sample_out_dir, threads, int(mem_mb/1024)], 
            cores=threads,
            mem=mem_mb,
            partition=partition,
            time=10080, # 7天
            name=f"assembly_{sample}"
        )
        
        assembly_info.append((sample, contigs))
        
    return assembly_info


def run_binning_workflow(workflow, assembly_info, clean_files_map, output_base_dir, 
                         threads=16, mem_mb=32000, partition=None, env_bowtie="", env_binning=""):
    """
    [Step 3 & 4] 整合流程: Build Index -> Mapping -> MetaBAT2 -> CheckM
    """
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
        index_target = index_base + ".1.bt2"
        sorted_bam = os.path.join(assembly_dir, f"{sample}.sorted.bam")
        
        log_build = os.path.join(assembly_dir, "bowtie2_build.log")

        # 2.1 Bowtie2 Build Index
        workflow.add_task_gridable(
            f"{env_bowtie} bowtie2-build --threads [args[0]] [depends[0]] [args[1]]",
            depends=contigs,
            targets=index_target,
            args=[threads, index_base],
            cores=threads,
            mem=mem_mb, 
            partition=partition,
            time=1440, # 24小时
            name=f"index_{sample}"
        )

        # --- 2.2 Bowtie2 Mapping ---
        # 管道命令的重定向：(cmd1 | cmd2) > log 2>&1
        log_map = os.path.join(assembly_dir, "bowtie2_mapping.log")
        map_mem = max(mem_mb, 32000) 
        cmd = (
            f"{env_bowtie} /bin/bash -c '"
            f"set -e; "
            f"bowtie2 -p [args[0]] -x [args[1]] -1 [depends[0]] -2 [depends[1]] | "
            f"samtools sort -@ [args[0]] -o [targets[0]] - && "
            f"samtools index [targets[0]]"
            f"' "
            f"> {log_map} 2>&1"
        )

        workflow.add_task_gridable(
            cmd,
            depends=[r1, r2, index_target],
            targets=sorted_bam,
            args=[threads, index_base],
            cores=threads,
            mem=map_mem,
            partition=partition,
            time=2880,
            name=f"map_{sample}"
        )

        # 3. MetaBAT2
        sample_bin_dir = os.path.join(dir_binning, sample)
        depth_file = os.path.join(sample_bin_dir, f"{sample}.depth.txt")
        metabat_done = os.path.join(sample_bin_dir, "metabat.done")
        log_metabat = os.path.join(sample_bin_dir, "metabat.log")
        
        workflow.add_task_gridable(
            f"mkdir -p [args[0]] && "
            f"{env_binning} /bin/bash -c '"
            f"set -e; "
            f"jgi_summarize_bam_contig_depths --outputDepth [targets[0]] [depends[0]] && "
            f"metabat2 -i [depends[1]] -a [targets[0]] -o [args[0]]/bin -t [args[1]] --minContig 1500"
            f"' "
            f"> {log_metabat} 2>&1 && "
            f"touch [targets[1]]",
            depends=[sorted_bam, contigs],
            targets=[depth_file, metabat_done],
            args=[sample_bin_dir, threads],
            cores=threads,
            mem=mem_mb,
            partition=partition,
            time=1440,
            name=f"binning_{sample}"
        )
        
        # 4. CheckM
        sample_checkm_dir = os.path.join(dir_checkm, sample)
        checkm_stats = os.path.join(sample_checkm_dir, "storage", "bin_stats_ext.tsv")
        checkm_mem = max(mem_mb, 64000)
        log_checkm = os.path.join(sample_checkm_dir, "checkm.log")
        
        workflow.add_task_gridable(
            f"mkdir -p [args[2]] && "
            f"{env_binning} checkm lineage_wf -t [args[0]] -x fa [args[1]] [args[2]] "
            f"> {log_checkm} 2>&1",
            depends=metabat_done,
            targets=checkm_stats,
            args=[threads, sample_bin_dir, sample_checkm_dir],
            cores=threads,
            mem=checkm_mem,
            partition=partition,
            time=2880, # 2天
            name=f"checkm_{sample}"
        )
        
        all_checkm_stats.append(checkm_stats)
        
    return all_checkm_stats

def _py_prepare_drep(task):
    """
    [Internal Helper] dRep 输入数据预处理
    """
    bins_link_dir = task.args[0]
    binning_base_dir = task.args[1]
    checkm_csv_path = task.targets[0].name
    
    if not os.path.exists(bins_link_dir):
        os.makedirs(bins_link_dir)
        
    all_parsed_data = []
    
    for checkm_file in task.depends:
        f_path = checkm_file.name
        try:
            sample_id = f_path.split(os.sep)[-3]
            df = pd.read_csv(f_path, sep='\t', header=None, names=['bin_id', 'dict_string'])
            for _, row in df.iterrows():
                stats = ast.literal_eval(row['dict_string'])
                bin_filename = f"{row['bin_id']}.fa"
                original_bin_path = os.path.join(binning_base_dir, sample_id, bin_filename)
                
                new_link_name = f"{sample_id}_{bin_filename}"
                link_path = os.path.join(bins_link_dir, new_link_name)
                
                if os.path.exists(original_bin_path):
                    if not os.path.exists(link_path):
                        os.symlink(os.path.abspath(original_bin_path), link_path)
                    
                    all_parsed_data.append({
                        'genome': new_link_name,
                        'completeness': stats['Completeness'],
                        'contamination': stats['Contamination']
                    })
        except Exception as e:
            print(f"Error parsing CheckM file {f_path}: {e}")

    if all_parsed_data:
        pd.DataFrame(all_parsed_data).to_csv(checkm_csv_path, index=False)
    else:
        print("Warning: No valid bins found for dRep.")
        open(checkm_csv_path, 'w').close()


def run_drep(workflow, checkm_results, drep_work_dir, binning_base_dir, 
             threads=32, mem_mb=200000, partition=None, env=""):
    """
    [Step 5] dRep 去冗余
    """
    bins_link_dir = os.path.join(drep_work_dir, "all_bins_linked")
    drep_info_csv = os.path.join(drep_work_dir, "checkm2drep.csv")
    
    workflow.add_task(
        _py_prepare_drep,
        depends=checkm_results,
        targets=drep_info_csv,
        args=[bins_link_dir, binning_base_dir],
        name="prepare_drep_inputs"
    )
    
    drep_out_genomes = os.path.join(drep_work_dir, "dereplicated_genomes")
    log_drep = os.path.join(drep_work_dir, "drep.log")
    
    workflow.add_task_gridable(
        f"{env} dRep dereplicate [args[0]] -p [args[1]] -g [args[2]]/*.fa "
        f"--genomeInfo [depends[0]] --completeness 50 --contamination 5 "
        f"-sa 0.95 -nc 0.3 -pa 0.90 "
        f"> {log_drep} 2>&1",
        depends=drep_info_csv,
        targets=TrackedDirectory(drep_out_genomes),
        args=[drep_work_dir, threads, bins_link_dir],
        cores=threads,
        mem=mem_mb,
        partition=partition,
        time=10080,
        name="drep_run"
    )
    
    return drep_out_genomes


def run_quantification(workflow, sgb_dir, clean_files_map, output_dir, 
                       threads=16, mem_mb=64000, partition=None, env_bowtie="", env_coverm=""):
    """
    [Step 6] 丰度计算
    """
    db_dir = os.path.join(output_dir, "database")
    sgb_fasta = os.path.join(db_dir, "SGBs.fna")
    sgb_index_prefix = os.path.join(db_dir, "SGBs_index")
    sgb_index_target = sgb_index_prefix + ".1.bt2"
    log_index = os.path.join(db_dir, "bowtie2_build.log")
    
    workflow.add_task_gridable(
        f"mkdir -p [args[0]] && "
        f"cat [depends[0]]/*.fa > [targets[0]] && "
        f"{env_bowtie} bowtie2-build --threads [args[1]] [targets[0]] [args[2]] "
        f"> {log_index} 2>&1",
        depends=TrackedDirectory(sgb_dir),
        targets=[sgb_fasta, sgb_index_target],
        args=[db_dir, threads, sgb_index_prefix],
        cores=threads,
        mem=mem_mb,
        partition=partition,
        time=1440,
        name="index_sgbs"
    )
    
    sgb_bams = []
    bam_dir = os.path.join(output_dir, "BAMs")
    
    for sample, (r1, r2) in clean_files_map.items():
        sgb_bam = os.path.join(bam_dir, f"{sample}.sorted.bam")
        log_map = os.path.join(bam_dir, f"{sample}_mapping.log")
        
        cmd = (
            f"mkdir -p [args[0]] && "
            f"{env_bowtie} /bin/bash -c '"
            f"set -e; "
            f"bowtie2 -p [args[1]] -x [args[2]] -1 [depends[0]] -2 [depends[1]] | "
            f"samtools sort -@ [args[1]] -o [targets[0]] - && "
            f"samtools index [targets[0]]"
            f"' "
            f"> {log_map} 2>&1"
        )
        
        workflow.add_task_gridable(
            cmd,
            depends=[r1, r2, sgb_index_target],
            targets=sgb_bam,
            args=[bam_dir, threads, sgb_index_prefix],
            cores=threads,
            mem=mem_mb,
            partition=partition,
            time=2880,
            name=f"map_sgb_{sample}"
        )
        sgb_bams.append(sgb_bam)
        
    abundance_table = os.path.join(output_dir, "sgb_abundance_matrix.tsv")
    log_coverm = os.path.join(output_dir, "coverm.log")
    
    workflow.add_task_gridable(
        f"{env_coverm} coverm genome -m relative_abundance --bam-files [args[0]]/*.sorted.bam "
        f"-f [depends[0]]/*.fa -t [args[1]] > [targets[0]] "
        f"2> {log_coverm}", # CoverM 结果在 stdout, 日志在 stderr
        depends=[TrackedDirectory(sgb_dir)] + sgb_bams, 
        targets=abundance_table,
        args=[bam_dir, threads],
        cores=threads,
        mem=mem_mb, 
        partition=partition,
        time=1440,
        name="coverm_calc"
    )
    
    return abundance_table


def run_gtdbtk(workflow, sgb_dir, output_dir, threads=32, mem_mb=500000, partition=None, env=""):
    """
    [Step 7] GTDB-Tk 物种注释
    """
    gtdb_summary = os.path.join(output_dir, "gtdbtk.bac120.summary.tsv")
    msa_file = os.path.join(output_dir, "align", "gtdbtk.bac120.user_msa.fasta")
    log_classify = os.path.join(output_dir, "gtdbtk_classify.log")
    
    # 1. Classify
    workflow.add_task_gridable(
        f"{env} gtdbtk classify_wf --genome_dir [depends[0]] --out_dir [args[0]] "
        f"--extension fa --cpus [args[1]] --skip_ani_screen "
        f"> {log_classify} 2>&1",
        depends=TrackedDirectory(sgb_dir),
        targets=[gtdb_summary, msa_file], 
        args=[output_dir, threads],
        cores=threads,
        mem=mem_mb,
        partition=partition,
        time=10080, 
        name="gtdbtk_classify"
    )
    
    # 2. Infer Tree
    tree_out_dir = os.path.join(output_dir, "classify")
    log_infer = os.path.join(output_dir, "gtdbtk_infer.log")
    
    workflow.add_task_gridable(
        f"{env} gtdbtk infer --out_dir [args[0]] --cpus [args[1]] --msa_file [depends[0]] "
        f"> {log_infer} 2>&1",
        depends=msa_file, 
        targets=TrackedDirectory(tree_out_dir),
        args=[tree_out_dir, threads],
        cores=threads,
        mem=128000, 
        partition=partition,
        time=2880,
        name="gtdbtk_infer"
    )