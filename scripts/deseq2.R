#!/usr/bin/env Rscript
# Differential expression with DESeq2.
#
# Produces, per contrast: a full results table (all genes, so nothing is hidden by an
# arbitrary threshold), and shrunken log2 fold changes for ranking/plotting. Also saves
# the DESeqDataSet and a variance-stabilised matrix for the report's PCA and heatmaps.

suppressPackageStartupMessages({
    library(DESeq2)
})

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output"); sink(log_con, type = "message")
on.exit({ sink(type = "message"); sink(type = "output"); close(log_con) }, add = TRUE)

counts_f   <- snakemake@input[["counts"]]
coldata_f  <- snakemake@input[["coldata"]]
design_str <- snakemake@params[["design"]]
contrasts  <- snakemake@params[["contrasts"]]
alpha      <- as.numeric(snakemake@params[["alpha"]])
lfc_thresh <- as.numeric(snakemake@params[["lfc_threshold"]])
min_count  <- as.numeric(snakemake@params[["min_count"]])
outdir     <- snakemake@params[["outdir"]]
# Optional ordered-dose trend test; empty list when not configured.
trend_cfg  <- tryCatch(snakemake@params[["dose_trend"]], error = function(e) list())
if (is.null(trend_cfg)) trend_cfg <- list()

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

counts <- read.delim(counts_f, row.names = 1, check.names = FALSE)
coldata <- read.delim(coldata_f, row.names = 1, check.names = FALSE,
                      stringsAsFactors = TRUE)

# DESeq2 requires the column order of the count matrix to match the row order of coldata
common <- intersect(colnames(counts), rownames(coldata))
if (length(common) != ncol(counts)) {
    stop("Sample mismatch between count matrix and sample sheet.\n",
         "  in counts only : ", paste(setdiff(colnames(counts), rownames(coldata)), collapse = ", "), "\n",
         "  in sheet only  : ", paste(setdiff(rownames(coldata), colnames(counts)), collapse = ", "))
}
counts <- counts[, common, drop = FALSE]
coldata <- coldata[common, , drop = FALSE]

# tximport's lengthScaledTPM counts are non-integer; DESeq2 wants integers
counts <- round(as.matrix(counts))
mode(counts) <- "integer"

message("Design: ~", design_str)
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = as.formula(paste("~", design_str)))

# Light pre-filter: drop genes with almost no signal. This mainly speeds things up and
# improves the multiple-testing correction; DESeq2's independent filtering does the
# statistically meaningful version of this later.
keep <- rowSums(counts(dds)) >= min_count
message(sprintf("Pre-filter: keeping %d / %d genes with total count >= %g",
                sum(keep), nrow(dds), min_count))
dds <- dds[keep, ]

dds <- DESeq(dds)
saveRDS(dds, snakemake@output[["dds"]])

# Variance-stabilised values for PCA / clustering / heatmaps. vst() needs enough genes
# to fit its dispersion trend; fall back to the slower-but-robust varianceStabilizingTransformation.
vsd <- tryCatch(
    vst(dds, blind = FALSE),
    error = function(e) {
        message("vst() failed (", conditionMessage(e), "); falling back to varianceStabilizingTransformation()")
        varianceStabilizingTransformation(dds, blind = FALSE)
    }
)
vst_mat <- assay(vsd)
write.table(data.frame(gene_id = rownames(vst_mat), vst_mat, check.names = FALSE),
            snakemake@output[["vst"]], sep = "\t", quote = FALSE, row.names = FALSE)

norm_counts <- counts(dds, normalized = TRUE)
write.table(data.frame(gene_id = rownames(norm_counts), round(norm_counts, 3),
                       check.names = FALSE),
            snakemake@output[["normalized"]], sep = "\t", quote = FALSE, row.names = FALSE)

# ------------------------------------------------------------------ contrasts ------
gene_names <- NULL
if (!is.null(snakemake@input[["tx2gene"]]) && length(snakemake@input[["tx2gene"]]) > 0) {
    t2g <- read.delim(snakemake@input[["tx2gene"]], stringsAsFactors = FALSE)
    if (all(c("gene_id", "gene_name") %in% colnames(t2g))) {
        gene_names <- unique(t2g[, c("gene_id", "gene_name")])
        gene_names <- gene_names[!duplicated(gene_names$gene_id), ]
        rownames(gene_names) <- gene_names$gene_id
    }
}

write_results <- function(res, name) {
    d <- as.data.frame(res)
    d$gene_id <- rownames(d)
    if (!is.null(gene_names)) {
        d$gene_name <- gene_names[d$gene_id, "gene_name"]
    } else {
        d$gene_name <- d$gene_id
    }
    cols <- c("gene_id", "gene_name", "baseMean", "log2FoldChange", "lfcSE",
              "stat", "pvalue", "padj")
    d <- d[, intersect(cols, colnames(d)), drop = FALSE]
    d <- d[order(d$padj, -abs(d$log2FoldChange), na.last = TRUE), ]

    dir.create(file.path(outdir, name), recursive = TRUE, showWarnings = FALSE)
    out_f <- file.path(outdir, name, "results.tsv")
    write.table(d, out_f, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

    sig <- subset(d, !is.na(padj) & padj < alpha & abs(log2FoldChange) >= lfc_thresh)
    write.table(sig, file.path(outdir, name, "significant.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    message(sprintf("%s: %d genes tested, %d significant (padj < %g, |log2FC| >= %g)",
                    name, nrow(d), nrow(sig), alpha, lfc_thresh))
    invisible(d)
}

if (length(contrasts) > 0) {
    # Each element arrives as "factor|numerator|denominator" (see analysis.smk)
    for (ct_str in as.character(contrasts)) {
        ct <- strsplit(ct_str, "|", fixed = TRUE)[[1]]
        if (length(ct) != 3) {
            stop("Malformed contrast ", sQuote(ct_str),
                 ": expected factor|numerator|denominator")
        }
        factor_name <- ct[1]; num <- ct[2]; den <- ct[3]
        name <- paste(factor_name, num, "vs", den, sep = "_")
        message("Contrast: ", name)
        res <- results(dds, contrast = c(factor_name, num, den), alpha = alpha)
        # apeglm needs a coefficient rather than a contrast; ashr works directly on
        # contrasts and is the appropriate choice here.
        res_shrunk <- tryCatch(
            lfcShrink(dds, contrast = c(factor_name, num, den), res = res, type = "ashr"),
            error = function(e) {
                message("lfcShrink failed (", conditionMessage(e), "); using unshrunken LFCs")
                res
            }
        )
        write_results(res_shrunk, name)
    }
} else {
    # No explicit contrast: report the last coefficient, which is DESeq2's default
    rn <- resultsNames(dds)
    coef <- tail(rn, 1)
    message("No contrasts configured; reporting coefficient: ", coef)
    res <- results(dds, name = coef, alpha = alpha)
    write_results(res, gsub("[^A-Za-z0-9._-]", "_", coef))
}

writeLines(resultsNames(dds), snakemake@output[["coefficients"]])
message("DESeq2 finished.")

# ------------------------------------------------------------- dose trend ---------
# A dose series is ordered information, and two pairwise contrasts throw that ordering
# away. Testing a single trend coefficient across the whole series asks the question the
# experiment was designed to answer ("does expression move with dose?") in one test
# rather than several, which means more power and one multiple-testing burden instead of
# one per comparison.
#
# Two encodings, because the right one depends on the spacing of the doses:
#
#   ordered  the dose levels are ranked and the rank is tested as a numeric covariate.
#            This is the linear trend an ordered factor's .L contrast would test, but
#            DESeq2 rejects ordered factors in a design formula outright ("the internal
#            steps do not work on ordered factors"), so the rank encoding is how that
#            test is actually obtained. Makes no assumption about the numeric spacing of
#            the levels, only their order. The safe default, and the right choice for a
#            vehicle / low / high design where the doses are not evenly spaced.
#
#   numeric  the dose column is coerced to a number and tested as a continuous
#            covariate, optionally on a log scale. More powerful when the doses really
#            are spaced the way the numbers say and the response is log-linear, which is
#            the usual pharmacological expectation. log1p rather than log because a
#            vehicle arm is dose zero.
#
# Either way this is a test for a monotone trend. It will miss a genuinely non-monotone
# response (a peak at the middle dose), which is why it is reported alongside the
# pairwise contrasts rather than instead of them.

if (length(trend_cfg) > 0 && !is.null(trend_cfg[["column"]])) {
    tcol <- as.character(trend_cfg[["column"]])
    encoding <- as.character(trend_cfg[["encoding"]] %||% "ordered")
    transform <- as.character(trend_cfg[["transform"]] %||% "log1p")

    if (!tcol %in% colnames(coldata)) {
        warning("dose_trend column ", sQuote(tcol), " is not in the sample sheet; ",
                "trend test skipped. Available: ", paste(colnames(coldata), collapse = ", "))
    } else {
        cd <- coldata
        raw <- as.character(cd[[tcol]])
        num <- suppressWarnings(as.numeric(raw))

        if (encoding == "numeric") {
            if (any(is.na(num))) {
                stop("dose_trend encoding 'numeric' needs column ", sQuote(tcol),
                     " to be numeric, but these values are not: ",
                     paste(unique(raw[is.na(num)]), collapse = ", "))
            }
            cd$.dose <- switch(transform,
                               none  = num,
                               log1p = log1p(num),
                               log10 = log10(num + 1),
                               stop("dose_trend transform must be none, log1p or log10"))
            coef_name <- ".dose"
        } else {
            # Order by numeric value where possible, otherwise by first appearance.
            lv <- unique(raw)
            lv <- if (any(is.na(num))) lv else lv[order(unique(num))]
            if (length(lv) < 3) {
                warning("dose_trend with only ", length(lv),
                        " levels is just a pairwise contrast; trend test skipped.")
                cd <- NULL
            } else {
                cd$.dose <- match(raw, lv)     # rank: 1, 2, 3, ... in dose order
            }
            coef_name <- ".dose"
        }

        if (!is.null(cd)) {
            # Keep any covariates from the configured design, dropping the dose factor
            # itself so it is not fitted twice under two encodings.
            other <- setdiff(trimws(strsplit(design_str, "\\+")[[1]]), c(tcol, ""))
            rhs <- paste(c(other, ".dose"), collapse = " + ")
            message("Dose trend: ~", rhs, "   (encoding: ", encoding, ")")

            dds_t <- DESeqDataSetFromMatrix(countData = counts, colData = cd,
                                            design = as.formula(paste("~", rhs)))
            dds_t <- dds_t[rowSums(counts(dds_t)) >= min_count, ]
            dds_t <- DESeq(dds_t, quiet = TRUE)

            nm <- resultsNames(dds_t)
            hit <- if (coef_name %in% nm) coef_name else grep("\\.dose", nm, value = TRUE)[1]
            if (is.na(hit)) {
                warning("could not find a dose coefficient among: ",
                        paste(nm, collapse = ", "))
            } else {
                res_t <- results(dds_t, name = hit, alpha = alpha)
                res_t <- lfcShrink(dds_t, coef = hit, res = res_t, type = "ashr")
                write_results(res_t, "dose_trend")
                message("Dose trend: ", sum(res_t$padj < alpha, na.rm = TRUE),
                        " genes at FDR < ", alpha, " on coefficient ", hit)
            }
        }
    }
}
