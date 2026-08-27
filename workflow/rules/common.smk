"""Shared setup: config defaults, sample sheet loading, pre-flight validation.

Everything that can be checked cheaply is checked here, at DAG-construction time, so
that misconfiguration fails in seconds rather than after hours of quantification.
"""

import os
import sys
from os.path import join

# Make scripts/ importable regardless of where snakemake was launched from.
# workflow.basedir is the directory holding the top-level Snakefile (i.e. workflow/),
# so the repo root is one level up. str() because Snakemake wraps this in its own
# source-file object rather than returning a plain path.
WORKFLOW_DIR = str(workflow.basedir)
REPO_ROOT = os.path.dirname(WORKFLOW_DIR)
SCRIPTS_DIR = join(REPO_ROOT, "scripts")
sys.path.insert(0, SCRIPTS_DIR)

from samplesheet import (  # noqa: E402
    load_samplesheet, validate_design, summarise, SampleSheetError,
)

# ---------------------------------------------------------------- config defaults ---
config.setdefault("outdir", "results")
config.setdefault("quantifier", "salmon")          # salmon | star_salmon
config.setdefault("trimmer", "fastp")              # fastp | none
config.setdefault("design", "condition")
config.setdefault("contrasts", [])
config.setdefault("reference", {})
config.setdefault("alpha", 0.05)                   # adjusted-p threshold
config.setdefault("lfc_threshold", 0.0)
config.setdefault("min_count", 10)                 # pre-filter: min total count per gene
config.setdefault("make_report", True)
config.setdefault("report_pdf", True)
config.setdefault("project_name", "Bulk RNA-seq analysis")
config.setdefault("report_top_n", 40)
config.setdefault("gsea", False)
config.setdefault("threads", {})
config["threads"].setdefault("trim", 4)
config["threads"].setdefault("index", 8)
config["threads"].setdefault("quant", 8)
config["threads"].setdefault("align", 8)

OUTDIR = config["outdir"]
QUANTIFIER = config["quantifier"]

_VALID_QUANTIFIERS = ("salmon", "star_salmon")
if QUANTIFIER not in _VALID_QUANTIFIERS:
    raise WorkflowError(
        f"config 'quantifier' must be one of {_VALID_QUANTIFIERS}, got {QUANTIFIER!r}"
    )

# ------------------------------------------------------------------- sample sheet ---
if "samplesheet" not in config:
    raise WorkflowError(
        "config is missing 'samplesheet'. Point it at a CSV/TSV with at least "
        "'sample' and 'fastq_1' columns - see config/samples_example.csv."
    )

try:
    SAMPLES, EXTRA_COLS = load_samplesheet(config["samplesheet"], check_files=True)
except SampleSheetError as e:
    raise WorkflowError(str(e))

SAMPLE_NAMES = list(SAMPLES)
if not SAMPLE_NAMES:
    raise WorkflowError(f"No samples found in {config['samplesheet']}")

# ------------------------------------------------------------ design / contrasts ---
# Differential expression is optional: with no usable design (e.g. no condition column)
# the pipeline still runs QC + quantification and reports those.
DE_PROBLEMS = validate_design(SAMPLES, EXTRA_COLS, config["design"],
                              config["contrasts"][0] if config["contrasts"] else None)
RUN_DE = not DE_PROBLEMS

for c in config["contrasts"][1:]:
    DE_PROBLEMS += validate_design(SAMPLES, EXTRA_COLS, config["design"], c)
RUN_DE = not DE_PROBLEMS

if config["contrasts"]:
    CONTRAST_NAMES = [f"{c[0]}_{c[1]}_vs_{c[2]}" for c in config["contrasts"]]
else:
    CONTRAST_NAMES = []

onstart:
    print("=" * 72)
    print(f"  {config['project_name']}")
    print("=" * 72)
    print(summarise(SAMPLES, EXTRA_COLS))
    print(f"  quantifier : {QUANTIFIER}")
    print(f"  output dir : {OUTDIR}")
    if RUN_DE:
        print(f"  design     : ~{config['design']}")
        print(f"  contrasts  : {CONTRAST_NAMES or '(all coefficients)'}")
    else:
        print("  differential expression: SKIPPED")
        for p in DE_PROBLEMS:
            print(f"    - {p}")
        print("    (QC and quantification will still run)")
    print("=" * 72)

# ------------------------------------------------------------------ path helpers ---
def sample_units(wildcards):
    return SAMPLES[wildcards.sample]["units"]


def is_paired(sample):
    return SAMPLES[sample]["paired"]


ALL_PAIRED = all(is_paired(s) for s in SAMPLE_NAMES)
ANY_PAIRED = any(is_paired(s) for s in SAMPLE_NAMES)


def raw_fastqs(wildcards, mate=1):
    """All FASTQ files for one sample and mate, in sample sheet order."""
    key = "fastq_1" if mate == 1 else "fastq_2"
    return [u[key] for u in SAMPLES[wildcards.sample]["units"] if u[key]]


def merged_fastq(sample, mate):
    return join(OUTDIR, "01_reads", f"{sample}_R{mate}.fastq.gz")


def trimmed_fastq(sample, mate):
    if config["trimmer"] == "none":
        return merged_fastq(sample, mate)
    return join(OUTDIR, "02_trimmed", f"{sample}_R{mate}.trimmed.fastq.gz")


def quant_input(wildcards):
    """FASTQs handed to the quantifier for one sample."""
    s = wildcards.sample
    if is_paired(s):
        return [trimmed_fastq(s, 1), trimmed_fastq(s, 2)]
    return [trimmed_fastq(s, 1)]


def _require_reference(*keys):
    ref = config["reference"]
    missing = [k for k in keys if not ref.get(k)]
    if missing:
        raise WorkflowError(
            f"quantifier '{QUANTIFIER}' needs reference file(s) {missing} set under "
            f"'reference:' in the config. See docs/usage.md and "
            f"scripts/download_references.sh."
        )
    for k in keys:
        p = ref[k]
        # index directories may be built by the pipeline; input files must exist now
        if k in ("transcriptome_fasta", "genome_fasta", "gtf") and not os.path.exists(p):
            raise WorkflowError(f"reference '{k}' does not exist: {p}")
    return [ref[k] for k in keys]


# Salmon index: use a prebuilt one if given, otherwise build from the transcriptome
SALMON_INDEX = config["reference"].get("salmon_index") or join(OUTDIR, "00_index", "salmon")
BUILD_SALMON_INDEX = not config["reference"].get("salmon_index")

STAR_INDEX = config["reference"].get("star_index") or join(OUTDIR, "00_index", "star")
BUILD_STAR_INDEX = not config["reference"].get("star_index")

if QUANTIFIER == "salmon":
    if BUILD_SALMON_INDEX:
        _require_reference("transcriptome_fasta")
    _require_reference("gtf")  # needed for tx2gene
elif QUANTIFIER == "star_salmon":
    if BUILD_STAR_INDEX:
        _require_reference("genome_fasta", "gtf")
    else:
        _require_reference("gtf")
    _require_reference("transcriptome_fasta")


# ----------------------------------------------------------------- target outputs ---
def final_outputs():
    out = [
        join(OUTDIR, "03_quant", "gene_counts.tsv"),
        join(OUTDIR, "03_quant", "gene_tpm.tsv"),
        join(OUTDIR, "03_quant", "transcript_counts.tsv"),
        join(OUTDIR, "00_qc", "multiqc_report.html"),
    ]
    if RUN_DE:
        out.append(join(OUTDIR, "04_de", "dds.rds"))
        out += expand(join(OUTDIR, "04_de", "{contrast}", "results.tsv"),
                      contrast=CONTRAST_NAMES) if CONTRAST_NAMES else []
    if config["make_report"]:
        out.append(join(OUTDIR, "05_report", "rnaseq_report.html"))
        if config["report_pdf"]:
            out.append(join(OUTDIR, "05_report", "rnaseq_report.pdf"))
    return out
