#!/usr/bin/env bash
# End-to-end smoke test on a small public dataset.
#
# Downloads ~28 MB of subsampled yeast RNA-seq from nf-core/test-datasets (GSE110004,
# a RAP1 induction experiment) plus a tiny reference, then runs the full pipeline:
# merge -> trim -> quantify -> tximport -> DESeq2 -> report.
#
# The sample sheet is deliberately awkward, because that is what makes it a useful test:
#   * WT_REP1 and RAP1_UNINDUCED_REP2 each span two sequencing units -> exercises merging
#   * the RAP1_UNINDUCED samples are single-end while the others are paired-end
#   * three conditions with two replicates each -> a real DESeq2 design
#
# Usage:
#   ./test/run_test.sh              # download if needed, then run
#   ./test/run_test.sh --dry-run    # build the DAG only, run nothing
#
# Runtime is a couple of minutes on a laptop; the whole dataset maps to ~124 genes.

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
REPO_ROOT="$PWD"

DRY_RUN=""
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN="-n"

# Snakemake manages the conda envs by default. Set USE_CONDA=0 if the tools are already
# on PATH (e.g. you created envs/environment.yml and envs/r.yml yourself).
CONDA_FLAG="--use-conda"
[[ "${USE_CONDA:-1}" == "0" ]] && CONDA_FLAG=""

DATA_DIR="test/data"
REF_DIR="$DATA_DIR/reference"
FQ_DIR="$DATA_DIR/fastq"
RAW="https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq"

mkdir -p "$REF_DIR" "$FQ_DIR"

fetch() {  # fetch <url> <dest>
    if [[ -s "$2" ]]; then return 0; fi
    echo "  downloading $(basename "$2")"
    curl -sfL -o "$2" "$1" || { echo "FAILED to download $1" >&2; exit 1; }
}

echo "==> Fetching reference (tiny yeast subset)"
fetch "$RAW/reference/genome.fa"            "$REF_DIR/genome.fa"
fetch "$RAW/reference/genes.gtf"            "$REF_DIR/genes.gtf"
fetch "$RAW/reference/transcriptome.fasta"  "$REF_DIR/transcriptome.fasta"

echo "==> Fetching FASTQ files"
for acc in SRR6357070_1 SRR6357070_2 SRR6357071_1 SRR6357071_2 SRR6357072_1 SRR6357072_2 \
           SRR6357073_1 SRR6357074_1 SRR6357075_1 SRR6357076_1 SRR6357076_2 \
           SRR6357077_1 SRR6357077_2; do
    fetch "$RAW/testdata/GSE110004/${acc}.fastq.gz" "$FQ_DIR/${acc}.fastq.gz"
done

echo "==> Writing sample sheet"
cat > "$DATA_DIR/samples.csv" <<EOF
sample,fastq_1,fastq_2,condition
WT_REP1,fastq/SRR6357070_1.fastq.gz,fastq/SRR6357070_2.fastq.gz,WT
WT_REP1,fastq/SRR6357071_1.fastq.gz,fastq/SRR6357071_2.fastq.gz,WT
WT_REP2,fastq/SRR6357072_1.fastq.gz,fastq/SRR6357072_2.fastq.gz,WT
RAP1_UNINDUCED_REP1,fastq/SRR6357073_1.fastq.gz,,RAP1_UNINDUCED
RAP1_UNINDUCED_REP2,fastq/SRR6357074_1.fastq.gz,,RAP1_UNINDUCED
RAP1_UNINDUCED_REP2,fastq/SRR6357075_1.fastq.gz,,RAP1_UNINDUCED
RAP1_IAA_30M_REP1,fastq/SRR6357076_1.fastq.gz,fastq/SRR6357076_2.fastq.gz,RAP1_IAA_30M
RAP1_IAA_30M_REP2,fastq/SRR6357077_1.fastq.gz,fastq/SRR6357077_2.fastq.gz,RAP1_IAA_30M
EOF

echo "==> Validating sample sheet before running anything"
python3 scripts/samplesheet.py "$DATA_DIR/samples.csv" condition condition,RAP1_IAA_30M,WT

echo "==> Writing test config"
cat > "$DATA_DIR/config_test.yaml" <<EOF
project_name: "Pipeline smoke test (yeast RAP1, subsampled)"
samplesheet: "$REPO_ROOT/$DATA_DIR/samples.csv"
outdir: "$REPO_ROOT/$DATA_DIR/results"
quantifier: "salmon"
trimmer: "fastp"
reference:
  transcriptome_fasta: "$REPO_ROOT/$REF_DIR/transcriptome.fasta"
  gtf: "$REPO_ROOT/$REF_DIR/genes.gtf"
design: "condition"
contrasts:
  - ["condition", "RAP1_IAA_30M", "WT"]
  - ["condition", "RAP1_UNINDUCED", "WT"]
alpha: 0.05
min_count: 5
make_report: true
report_pdf: false
threads:
  trim: 2
  index: 2
  quant: 2
EOF

echo "==> Running pipeline"
snakemake -s workflow/Snakefile \
    --configfile "$DATA_DIR/config_test.yaml" \
    --cores "${CORES:-4}" \
    $CONDA_FLAG \
    $DRY_RUN

if [[ -z "$DRY_RUN" ]]; then
    echo
    echo "==> Checking expected outputs"
    fail=0
    for f in "$DATA_DIR/results/03_quant/gene_counts.tsv" \
             "$DATA_DIR/results/04_de/dds.rds" \
             "$DATA_DIR/results/04_de/condition_RAP1_IAA_30M_vs_WT/results.tsv" \
             "$DATA_DIR/results/05_report/rnaseq_report.html"; do
        if [[ -s "$f" ]]; then echo "  OK   $f"; else echo "  MISS $f"; fail=1; fi
    done
    echo
    if [[ $fail -eq 0 ]]; then
        n=$(($(wc -l < "$DATA_DIR/results/03_quant/gene_counts.tsv") - 1))
        echo "Smoke test PASSED - quantified $n genes across 6 samples."
    else
        echo "Smoke test FAILED - see logs under $DATA_DIR/results/logs/" >&2
        exit 1
    fi
fi
