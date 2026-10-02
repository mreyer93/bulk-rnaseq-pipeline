#!/usr/bin/env Rscript
# Render the analysis report from a Snakemake rule.
#
# HTML and PDF are rendered by separate rules, so a broken LaTeX toolchain costs you the
# PDF but never the HTML. Each format renders entirely inside its own scratch directory
# (figures, LaTeX files and the report itself), and only the finished report is moved
# into place. Snakemake runs the two rules concurrently; rendered straight into the shared
# output directory, both wrote figures into the same <name>_files/ folder, which rmarkdown
# deletes when a render finishes.

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output"); sink(log_con, type = "message")
on.exit({ sink(type = "message"); sink(type = "output"); close(log_con) }, add = TRUE)

rmd_in   <- normalizePath(snakemake@params[["rmd"]], mustWork = TRUE)
out_file <- snakemake@output[[1]]
fmt      <- snakemake@params[["format"]]

out_dir <- dirname(out_file)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, mustWork = TRUE)

# Removed once the report is in place (kept after a failure, for debugging). Not via
# on.exit(): Snakemake runs this file at top level, where on.exit() never fires.
inter_dir <- file.path(out_dir, paste0(".render_", fmt))
unlink(inter_dir, recursive = TRUE)
dir.create(inter_dir, recursive = TRUE, showWarnings = FALSE)

p <- snakemake@params
report_params <- list(
    project_name  = p[["project_name"]],
    outdir        = normalizePath(p[["outdir"]], mustWork = FALSE),
    scripts_dir   = normalizePath(p[["scripts_dir"]], mustWork = TRUE),
    design        = p[["design"]],
    contrasts     = paste(p[["contrasts"]], collapse = ","),
    alpha         = p[["alpha"]],
    lfc_threshold = p[["lfc_threshold"]],
    top_n         = p[["top_n"]],
    quantifier    = p[["quantifier"]],
    run_de        = isTRUE(p[["run_de"]]),
    de_problems   = paste(p[["de_problems"]], collapse = " | ")
)

message("Rendering ", rmd_in, " as ", fmt)
rmarkdown::render(
    input             = rmd_in,
    output_format     = fmt,
    output_file       = basename(out_file),
    output_dir        = inter_dir,
    intermediates_dir = inter_dir,
    knit_root_dir     = getwd(),
    params            = report_params,
    envir             = new.env(parent = globalenv()),
    quiet             = FALSE
)

rendered <- file.path(inter_dir, basename(out_file))
if (!file.exists(rendered)) {
    stop("rmarkdown::render finished but ", rendered, " was not created")
}
if (!file.rename(rendered, out_file)) {
    stop("could not move ", rendered, " to ", out_file)
}
unlink(inter_dir, recursive = TRUE)
message("Wrote ", out_file)
