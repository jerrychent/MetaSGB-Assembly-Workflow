# MetaSGB Assembly Workflow

An anadama2 workflow for metagenomic SGB/MAG assembly and dereplication. The pipeline assumes all required tools are available from the currently active runtime environment. The KneadData database can be selected per run with `--kneaddata-db`.

## Workflow steps

1. `KneadData`: clean and decontaminate paired-end reads into `02_Cleandata/<sample>/`.
2. `MetaSPAdes`: assemble each sample into `03_Assembly/<sample>/contigs.fasta`.
3. `Bowtie2/Samtools`: map clean reads back to contigs and create sorted BAM files.
4. `MetaBAT2/CheckM`: bin contigs and estimate MAG quality.
5. `dRep`: combine all bins, use CheckM quality metadata, and dereplicate SGB representatives.
6. `CoverM`: map reads to dereplicated SGBs and produce an abundance matrix.
7. `GTDB-Tk`: classify SGB representatives and optionally infer a tree.

## Inputs

Raw paired-end reads are discovered by anadama2 from the input directory using `--extension-paired`. The default suffixes are:

```text
_1.fastq.gz,_2.fastq.gz
```

Sample IDs are inferred by removing the R1 suffix from each R1 filename.

## Example

```bash
python run_sgb.py \
  --input /path/to/raw_reads \
  --output /path/to/sgb_output \
  --resource-cfg config/resources.cfg \
  --extension-paired _1.fastq.gz,_2.fastq.gz \
  --kneaddata-db /path/to/kneaddata_db \
  --grid-jobs 40
```

The command does not expose per-tool environment prefixes. Activate the integrated environment before launching the workflow, and use `--kneaddata-db` to select the host/contaminant database.

## Bypass options

- `--bypass-kneaddata`: use paired clean reads under `02_Cleandata/**`. Both `*_paired_1.fastq` and `*_paired_1.fastq.gz` are supported.
- `--bypass-assembly`: use existing `03_Assembly/<sample>/contigs.fasta`.
- `--bypass-binning`: use existing `05_Checkm/<sample>/storage/bin_stats_ext.tsv`.
- `--bypass-drep`: use existing `06_dRep_SGBs/dereplicated_genomes/*.fa`.
- `--bypass-quant`: skip CoverM abundance calculation.
- `--bypass-gtdbtk`: skip GTDB-Tk classification and tree inference.

## Resource and tool parameters

All scheduler resources and tunable tool options are configured in `config/resources.cfg`.

Key sections:

- `[kneaddata]`: `threads`, `memory_mb`, `partition`, `sequencer_source`, FastQC toggles, and `extra_args`.
- `[assembly]`: MetaSPAdes resources and `extra_args`.
- `[binning]`: Bowtie2/MetaBAT2 resources, `min_contig`, `bowtie2_extra_args`, and `metabat_extra_args`.
- `[drep]`: dRep resources and thresholds: `completeness`, `contamination`, `secondary_ani`, `primary_ani`, `coverage`, and `extra_args`.
- `[quantification]`: CoverM resources, `coverm_method`, and `coverm_extra_args`.
- `[gtdbtk]`: GTDB-Tk resources, `marker_set`, `skip_ani_screen`, `infer_tree`, `infer_memory_mb`, and `extra_args`.

## Main outputs

- Clean reads: `02_Cleandata/<sample>/<sample>_paired_*.fastq`
- Assemblies: `03_Assembly/<sample>/contigs.fasta`
- Per-sample bins: `04_Binning/<sample>/bin.*.fa`
- CheckM stats: `05_Checkm/<sample>/storage/bin_stats_ext.tsv`
- Dereplicated SGBs: `06_dRep_SGBs/dereplicated_genomes/*.fa`
- Abundance matrix: `08_Abundance/sgb_abundance_matrix.tsv`
- GTDB-Tk summary: `07_GTDBTk/gtdbtk.<marker_set>.summary.tsv`

## Notes

- Bowtie2 indexes are tracked with `.done` sentinel files so both normal `.bt2` and large `.bt2l` indexes are supported.
- Mapping tasks use `set -euo pipefail` so Bowtie2 failures inside pipes propagate correctly.
- dRep input preparation raises a clear error if no valid bins are found.

