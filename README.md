# Bulk RNA-seq pipeline

A Snakemake pipeline taking bulk RNA-seq from raw FASTQ to differential expression and a
client-ready report. Tool choices follow [nf-core/rnaseq](https://nf-co.re/rnaseq), the
community reference implementation, so results are comparable with the standard:

```
fastp -> Salmon (or STAR + Salmon) -> tximport -> DESeq2 -> HTML/PDF report
```

## Two paths, same outputs

| | `quantifier: salmon` | `quantifier: star_salmon` |
|---|---|---|
| Method | Selective alignment from FASTQ | STAR genomic alignment, then Salmon |
| RAM (human) | ~4-8 GB | ~30-40 GB |
| Extra outputs | - | coordinate-sorted genomic BAM, splice junctions |
| Intended for | laptops, small servers | cloud VMs, compute servers |
| Config | [`config/config_local.yaml`](config/config_local.yaml) | [`config/config_cloud.yaml`](config/config_cloud.yaml) |

Both converge on the same gene counts, so downstream analysis, figures and the report are
identical. The split exists because a STAR human index does not fit in 16 GB of RAM —
this is the single biggest practical constraint in bulk RNA-seq on a laptop.

## Quickstart

```bash
mamba env create -f envs/environment.yml -n bulk-rnaseq
conda activate bulk-rnaseq

# check your sample sheet and design before running anything
python3 scripts/samplesheet.py my_samples.csv condition condition,treated,control

cp config/config_local.yaml my_config.yaml   # edit paths, design, contrasts
snakemake --configfile my_config.yaml --use-conda --cores 4 -n   # dry run
snakemake --configfile my_config.yaml --use-conda --cores 4
```

Verify the whole thing works first on a small public dataset (~28 MB, a couple of minutes):

```bash
./test/run_test.sh
```

## Sample sheets

CSV or TSV. The only hard requirements are a sample-name column and a first-FASTQ column;
common naming conventions are recognised automatically:

| Meaning | Accepted column names |
|---|---|
| Sample name | `sample`, `sample_id`, `Sample_ID`, `sample_name`, `name`, `id` |
| Read 1 | `fastq_1`, `fastq1`, `R1`, `read1`, `forward`, `fq1`, `fastq` |
| Read 2 | `fastq_2`, `fastq2`, `R2`, `read2`, `reverse`, `fq2` |
| Group | `condition`, `group`, `treatment`, `class`, `phenotype`, `status` |
| Batch | `batch`, `run`, `lane`, `flowcell`, `plate` |

Any other columns are carried through as metadata and can be used in the DESeq2 design.
Leave `fastq_2` blank for single-end. Repeat a sample name across rows to merge multiple
lanes or runs. Paths may be absolute or relative to the sample sheet.

`scripts/samplesheet.py` validates all of this — missing files, duplicated FASTQs across
samples, groups without replication, design terms that do not exist — and reports every
problem with a line number, before any compute is spent.

## Differential expression

Set a design and contrasts in the config:

```yaml
design: "batch + condition"          # variable of interest last
contrasts:
  - ["condition", "treated", "control"]   # log2FC is treated relative to control
```

If no usable design is present (no metadata columns, or only one group), the pipeline
still runs QC and quantification and says so in the report, rather than failing.

## Outputs

```
00_qc/     MultiQC report, fastp and FastQC output
03_quant/  gene_counts.tsv, gene_tpm.tsv, transcript_counts.tsv, per-sample Salmon output
04_de/     dds.rds, normalized_counts.tsv, vst.tsv, <contrast>/results.tsv + significant.tsv
05_report/ rnaseq_report.html / .pdf, and summary_tables/ backing every figure
logs/      one log per rule
```

See [docs/usage.md](docs/usage.md) for configuration detail and
[docs/outputs.md](docs/outputs.md) for what each file contains.

## Running in the cloud

[`cloud/gcp/`](cloud/gcp/README.md) provisions a VM sized for the `star_salmon` path and
a persistent data disk for references, with start/stop scripts to control cost.

## References

- [nf-core/rnaseq](https://nf-co.re/rnaseq) — the reference implementation these tool choices follow
- [Salmon](https://salmon.readthedocs.io/) — Patro et al., *Nat Methods* 2017
- [STAR](https://github.com/alexdobin/STAR) — Dobin et al., *Bioinformatics* 2013
- [tximport](https://bioconductor.org/packages/tximport/) — Soneson et al., *F1000Research* 2015
- [DESeq2](https://bioconductor.org/packages/DESeq2/) — Love et al., *Genome Biology* 2014
