"""Differential expression and the client-facing report."""

import csv as _csv


rule coldata:
    """Write the sample metadata table (the DESeq2 colData) from the sample sheet.

    Materialised as a file rather than passed in memory so that the DE step is
    reproducible standalone and the exact design inputs are recorded with the results.
    """
    output: join(OUTDIR, "04_de", "coldata.tsv")
    log: join(OUTDIR, "logs", "coldata.log")
    run:
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        cols = ["sample"] + [c for c in EXTRA_COLS]
        with open(output[0], "w", newline="") as fh:
            w = _csv.writer(fh, delimiter="\t")
            w.writerow(cols)
            for name, s in SAMPLES.items():
                w.writerow([name] + [s.get(c, "") for c in EXTRA_COLS])
        with open(log[0], "w") as fh:
            fh.write(f"Wrote colData for {len(SAMPLES)} samples: {cols}\n")


rule deseq2:
    input:
        counts = join(OUTDIR, "03_quant", "gene_counts.tsv"),
        coldata = join(OUTDIR, "04_de", "coldata.tsv"),
        tx2gene = join(OUTDIR, "00_index", "tx2gene.tsv"),
    output:
        dds = join(OUTDIR, "04_de", "dds.rds"),
        vst = join(OUTDIR, "04_de", "vst.tsv"),
        normalized = join(OUTDIR, "04_de", "normalized_counts.tsv"),
        coefficients = join(OUTDIR, "04_de", "coefficients.txt"),
        # The model is fitted once and every configured contrast is written in the same
        # pass, so all of them are declared here rather than via a separate rule.
        results = expand(join(OUTDIR, "04_de", "{contrast}", "results.tsv"),
                         contrast=CONTRAST_NAMES),
        significant = expand(join(OUTDIR, "04_de", "{contrast}", "significant.tsv"),
                             contrast=CONTRAST_NAMES),
    log: join(OUTDIR, "logs", "deseq2.log")
    params:
        design = config["design"],
        # Passed as "factor|numerator|denominator" strings rather than nested lists:
        # Snakemake flattens nested Python lists into a plain character vector when it
        # hands params to R, which silently destroys the grouping.
        contrasts = ["|".join(c) for c in config["contrasts"]],
        alpha = config["alpha"],
        lfc_threshold = config["lfc_threshold"],
        min_count = config["min_count"],
        outdir = join(OUTDIR, "04_de"),
    conda: "../../envs/r.yml"
    script: "../../scripts/deseq2.R"


def report_inputs(wildcards):
    deps = [
        join(OUTDIR, "03_quant", "gene_counts.tsv"),
        join(OUTDIR, "03_quant", "gene_tpm.tsv"),
        join(OUTDIR, "00_qc", "multiqc_report.html"),
    ]
    if RUN_DE:
        deps += [
            join(OUTDIR, "04_de", "dds.rds"),
            join(OUTDIR, "04_de", "vst.tsv"),
            join(OUTDIR, "04_de", "coldata.tsv"),
        ]
        deps += expand(join(OUTDIR, "04_de", "{contrast}", "results.tsv"),
                       contrast=CONTRAST_NAMES)
    return deps


_report_params = dict(
    scripts_dir = SCRIPTS_DIR,
    outdir = OUTDIR,
    project_name = config["project_name"],
    design = config["design"],
    contrasts = CONTRAST_NAMES,
    alpha = config["alpha"],
    lfc_threshold = config["lfc_threshold"],
    top_n = config["report_top_n"],
    quantifier = QUANTIFIER,
    run_de = RUN_DE,
    de_problems = DE_PROBLEMS,
)


rule report_html:
    input: report_inputs
    output: join(OUTDIR, "05_report", "rnaseq_report.html")
    log: join(OUTDIR, "logs", "report_html.log")
    params: rmd = join(SCRIPTS_DIR, "rnaseq_report.Rmd"), format = "html_document",
            **_report_params
    conda: "../../envs/r.yml"
    script: "../../scripts/render_report.R"


rule report_pdf:
    input: report_inputs
    output: join(OUTDIR, "05_report", "rnaseq_report.pdf")
    log: join(OUTDIR, "logs", "report_pdf.log")
    params: rmd = join(SCRIPTS_DIR, "rnaseq_report.Rmd"), format = "pdf_document",
            **_report_params
    conda: "../../envs/r.yml"
    script: "../../scripts/render_report.R"
