#!/usr/bin/env Rscript
# Render the analysis report from a Snakemake rule.
#
# HTML and PDF are rendered by separate rules using separate intermediates directories,
# so the two formats cannot clobber each other's temporary files and a broken LaTeX
# toolchain costs you the PDF but never the HTML.

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output"); sink(log_con, type = "message")
on.exit({ sink(type = "message"); sink(type = "output"); close(log_con) }, add = TRUE)

rmd_in   <- normalizePath(snakemake@params[["rmd"]], mustWork = TRUE)
out_file <- snakemake@output[[1]]
fmt      <- snakemake@params[["format"]]

out_dir <- dirname(out_file)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, mustWork = TRUE)

inter_dir <- file.path(out_dir, paste0(".render_", fmt))
dir.create(inter_dir, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(inter_dir, recursive = TRUE), add = TRUE)

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
    output_dir        = out_dir,
    intermediates_dir = inter_dir,
    knit_root_dir     = getwd(),
    params            = report_params,
    envir             = new.env(parent = globalenv()),
    quiet             = FALSE
)

if (!file.exists(out_file)) {
    stop("rmarkdown::render finished but ", out_file, " was not created")
}
message("Wrote ", out_file)
