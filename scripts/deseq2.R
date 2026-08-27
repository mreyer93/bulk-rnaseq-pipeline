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
