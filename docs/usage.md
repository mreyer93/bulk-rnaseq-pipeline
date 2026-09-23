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

## 2b. Transposable element annotation (optional)

Needed only if you want to quantify transposable elements and endogenous retroviruses
alongside genes. Skip this if you only care about gene-level expression.

`scripts/make_te_gtf.py` builds a TEtranscripts-compatible GTF from the UCSC RepeatMasker
track for any assembly UCSC hosts:

```
python scripts/make_te_gtf.py --genome mm39 --out references/mm39_TE.gtf
python scripts/make_te_gtf.py --genome hg38 --out references/hg38_TE.gtf
```

It is generated rather than downloaded on purpose. The prebuilt TE GTFs that TEtranscripts
documents are distributed by hand, and the lab file-share URL most tutorials cite now
returns 404. UCSC's RepeatMasker track is a stable, versioned source, so building from it
keeps the annotation reproducible.

By default it keeps LINE, SINE, LTR, DNA, Retroposon and RC elements, and drops simple
repeats, low-complexity regions, satellites and the small-RNA classes. Those are
repetitive but are not transposable elements, and including them inflates the count matrix
with features nobody will interpret. Override with `--classes` (or `--classes all`).

Check the build without running the pipeline:

```
python scripts/make_te_gtf.py --selftest
```

On mouse mm39 the output should contain the ERV families that matter for
chromatin-repression work, including the IAP elements (`IAPEz-int`, `IAPLTR*`), `MusD`/`ETn`
and `MMERVK10C-int`.

*One alignment caveat that decides whether TE analysis is even possible.* TE quantification
depends on multi-mapping reads, because young high-copy elements are near-identical across
loci. TEtranscripts' authors recommend `--winAnchorMultimapNmax 100 --outFilterMultimapNmax 100`
for STAR. STAR's default is 10, so BAMs produced by a standard pipeline (including most
sequencing-core deliverables) have already discarded the reads a TE analysis needs. If you
are handed BAMs rather than FASTQs, check the STAR command in the BAM header before
promising a TE result.

## 2c. Starting from a vendor's QC report

If a sequencing provider handed over MultiQC HTML and little else, the per-sample numbers
are still recoverable. MultiQC 1.x compresses its plot data inside the page, so the values
cannot be grepped and the report looks like a dead end. It is not:

```
python scripts/multiqc_extract.py --list multiqc_report.html
python scripts/multiqc_extract.py multiqc_report.html --outdir metrics/
```

One TSV per plot, samples as rows and series as columns, with a `total` column added. This
is enough to audit alignment rates, library sizes and featureCounts assignment across a
study before any of the primary data arrives, and it is often enough to settle
strandedness: run the provider's counts three ways and the correct setting assigns an
order of magnitude more reads than the wrong one.

Needs `pip install lzstring` for the compressed flavour.

## 2d. Three analyses beyond a standard contrast

All optional and all off by default.

*Transposable elements and ERVs.* Set `te_analysis.enabled: true` with a TE GTF from
section 2b. This aligns separately with STAR at `--outFilterMultimapNmax 100` and counts
genes and TEs together with TEcount, so both share one size factor. Outputs land in
`06_te/`: a combined matrix, gene-only and TE-only matrices, and a per-class/per-family
summary. `te_analysis.exclude_bed` drops regions from the BAM first, which matters more
than it sounds: a gene-free, repeat-dense window that is duplicated in the assembly will
otherwise be reported as a large LINE or LTR signal that is an artefact.

*Ordered dose trend.* Set `dose_trend: {column: dose, encoding: ordered}`. A dose series
carries ordering that two pairwise contrasts discard, and one trend coefficient across
the whole series is a single test rather than several. Note DESeq2 rejects ordered
factors in a design formula, so `ordered` ranks the levels and fits the rank as a numeric
covariate, which is the same linear trend an ordered factor's `.L` contrast would give.
Use `encoding: numeric` with `transform: log1p` when the doses really are spaced as the
numbers say.

Benchmarked rather than assumed: on a nine-dose, five-tissue mouse study the trend test
found more genes than the top-dose-versus-vehicle contrast in two tissues and fewer in
three, though one of those three was within 4% and is better read as a tie. Each test
found genes the other missed in every tissue.
Run both. Neither dominates, because they answer different questions: monotone movement
across the series versus a difference at the top dose.

*Cross-tissue concordance.* `scripts/concordance.R` finds genes moving the same way in
two independent sets of contrasts, which is the analysis behind a peripheral-biomarker
question ("what changes in blood when this changes in brain?").

```
Rscript scripts/concordance.R --a brain1.tsv,brain2.tsv --b blood1.tsv,blood2.tsv --outdir out/
```

It does not simply intersect two FDR-filtered lists. The intersection of two 5% FDR lists
has no controlled error rate, and reporting it as though it did is common. The primary
result is an intersection-union test: for a conjunction hypothesis the valid statistic is
the *maximum* of the component p-values, not a combination of them, so Fisher's or
Stouffer's method is the wrong tool here and will call genes that move in opposite
directions. Effect-size correlation and directional agreement are reported alongside,
because the hit count on its own does not measure shared biology: two tissues with
thousands of responding genes overlap by arithmetic.

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
