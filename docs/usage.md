# Usage

## 1. Install

```bash
mamba env create -f envs/environment.yml -n bulk-rnaseq
conda activate bulk-rnaseq
```

The R/Bioconductor stack (DESeq2, tximport, report rendering) lives in a second
environment, `envs/r.yml`. With `--use-conda` Snakemake creates and uses both
automatically; otherwise create it yourself and make sure `Rscript` is on `PATH`.

## 2. References

You need, from **the same annotation release**:

| Path in config | `salmon` | `star_salmon` | What it is |
|---|---|---|---|
| `reference.transcriptome_fasta` | required | required | cDNA/transcript FASTA |
| `reference.gtf` | required | required | gene annotation (transcript→gene mapping) |
| `reference.genome_fasta` | – | required | genome FASTA (for the STAR index) |

Mixing releases is the classic silent failure in RNA-seq: transcript IDs stop matching,
and you get an all-zero count matrix rather than an error. The pipeline guards against
this — `tximport_counts.R` reports the transcript-ID overlap and aborts if it is zero.

Ensembl, for example:

```bash
SP=homo_sapiens; REL=112
curl -O https://ftp.ensembl.org/pub/release-$REL/fasta/$SP/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
curl -O https://ftp.ensembl.org/pub/release-$REL/gtf/$SP/Homo_sapiens.GRCh38.$REL.gtf.gz
curl -O https://ftp.ensembl.org/pub/release-$REL/fasta/$SP/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```

If transcript IDs carry version suffixes in one file but not the other
(`ENST00000123.4` vs `ENST00000123`), set `strip_tx_version: true`.

## 3. Sample sheet

CSV or TSV; the separator is detected. Only a sample-name column and a first-FASTQ column
are mandatory. See the table in the [README](../README.md#sample-sheets) for accepted
column-name aliases.

```csv
sample,fastq_1,fastq_2,condition,batch
CTRL_1,fastq/C1_R1.fq.gz,fastq/C1_R2.fq.gz,control,b1
CTRL_2,fastq/C2_R1.fq.gz,fastq/C2_R2.fq.gz,control,b2
TREAT_1,fastq/T1_R1.fq.gz,fastq/T1_R2.fq.gz,treated,b1
TREAT_2,fastq/T2_R1.fq.gz,fastq/T2_R2.fq.gz,treated,b2
```

- **Single-end**: leave `fastq_2` empty. Single- and paired-end samples can coexist.
- **Multiple lanes/runs**: repeat the sample name on more than one row; they are concatenated.
- **Paths**: absolute, or relative to the sample sheet's own directory.
- **Extra columns** are carried through and usable in the design.

Always validate before running:

```bash
python3 scripts/samplesheet.py samples.csv "batch + condition" condition,treated,control
```

This checks that every FASTQ exists, that no FASTQ is shared between samples, that no
group lacks replication, and that every design term and contrast level actually exists.

## 4. Configure

Copy `config/config_local.yaml` (Salmon, laptop-sized) or `config/config_cloud.yaml`
(STAR+Salmon, server-sized) and edit. Key options:

| Option | Meaning |
|---|---|
| `quantifier` | `salmon` (~4-8 GB RAM) or `star_salmon` (~30-40 GB RAM) |
| `design` | DESeq2 formula over sample sheet columns; variable of interest **last** |
| `contrasts` | list of `[factor, numerator, denominator]`; log2FC is numerator vs denominator |
| `alpha` | adjusted-p threshold (default 0.05) |
| `lfc_threshold` | optional additional \|log2FC\| cut |
| `min_count` | drop genes below this total count before testing |
| `report_pdf` | render PDF as well as HTML (needs tectonic) |
| `threads` | per-stage thread counts; match your machine |

### Choosing a design

`design: "condition"` compares groups directly. `design: "batch + condition"` adjusts for
a batch effect while testing condition — put the variable of interest last, because
DESeq2's default coefficient is the last term. A batch that is perfectly confounded with
condition cannot be adjusted for by any model; the validator will not catch that, so
check your experimental design.

## 5. Run

```bash
snakemake --configfile my_config.yaml --use-conda --cores 8 -n   # dry run first
snakemake --configfile my_config.yaml --use-conda --cores 8
```

Useful flags: `-k` keep going past failures, `--rerun-incomplete` after an interrupted run,
`-p` print shell commands, `--report report.html` for Snakemake's own provenance report.

## Troubleshooting

| Symptom | Likely cause |
|---|---|
| tximport aborts on transcript-ID overlap | transcriptome FASTA and GTF from different releases; or set `strip_tx_version` |
| All-zero / near-zero counts | wrong reference species, or failed trimming — check the Salmon mapping rate in `logs/salmon_quant/` |
| STAR killed / out of memory | switch to `quantifier: salmon`, or use a machine with ≥32 GB RAM |
| "group(s) with no replication" | DESeq2 cannot estimate dispersion from n=1; add replicates or drop that group |
| Low mapping rate in one sample | contamination, adapter/quality issues, or degraded input; inspect its fastp and FastQC output |
| PDF report fails, HTML fine | LaTeX/tectonic missing — set `report_pdf: false`, or install `envs/r.yml` fully |

Every rule writes its own log under `<outdir>/logs/`. On failure, Snakemake names the
failing rule and its log path; that log almost always contains the actual cause.
