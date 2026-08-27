#!/usr/bin/env Rscript
# Aggregate per-sample Salmon quantifications to gene level with tximport.
#
# countsFromAbundance = "lengthScaledTPM" produces counts that are corrected for both
# library size and average transcript length, which is what makes them valid input to
# DESeq2 while remaining comparable between samples. See the tximport vignette.

suppressPackageStartupMessages({
    library(tximport)
})

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output"); sink(log_con, type = "message")
on.exit({ sink(type = "message"); sink(type = "output"); close(log_con) }, add = TRUE)

samples   <- snakemake@params[["samples"]]
quant_dir <- snakemake@params[["quant_dir"]]
tx2gene_f <- snakemake@input[["tx2gene"]]

files <- file.path(quant_dir, samples, "quant.sf")
names(files) <- samples
missing <- files[!file.exists(files)]
if (length(missing) > 0) {
    stop("Missing Salmon output for: ", paste(names(missing), collapse = ", "))
}

tx2gene <- read.delim(tx2gene_f, stringsAsFactors = FALSE)
if (!all(c("transcript_id", "gene_id") %in% colnames(tx2gene))) {
    stop("tx2gene table must have transcript_id and gene_id columns; got: ",
         paste(colnames(tx2gene), collapse = ", "))
}

# Sanity-check that the annotation and the index actually describe the same
# transcriptome. A near-total mismatch here is the classic silent failure in RNA-seq -
# it yields an all-zero count matrix rather than an error - so fail loudly instead.
first <- read.delim(files[1], stringsAsFactors = FALSE, nrows = 20000)
overlap <- length(intersect(first$Name, tx2gene$transcript_id))
frac <- overlap / min(nrow(first), nrow(tx2gene))
message(sprintf("Transcript ID overlap between quant.sf and tx2gene: %d (%.1f%%)",
                overlap, 100 * frac))
if (overlap == 0) {
    stop("No transcript IDs in common between the Salmon index and the GTF-derived ",
         "tx2gene table.\n",
         "  quant.sf examples : ", paste(head(first$Name, 3), collapse = ", "), "\n",
         "  tx2gene examples  : ", paste(head(tx2gene$transcript_id, 3), collapse = ", "), "\n",
         "The transcriptome FASTA and the GTF must come from the same annotation ",
         "release. If one carries version suffixes (ENST00000123.4) and the other ",
         "does not, set 'strip_tx_version: true' in the config.")
} else if (frac < 0.5) {
    warning("Fewer than half of transcripts match between the index and the ",
            "annotation - check they are the same release.")
}

txi <- tximport(files, type = "salmon", tx2gene = tx2gene[, c("transcript_id", "gene_id")],
                countsFromAbundance = "lengthScaledTPM",
                ignoreTxVersion = isTRUE(snakemake@config[["strip_tx_version"]]))

txi_tx <- tximport(files, type = "salmon", txOut = TRUE)

write_mat <- function(mat, path, idname) {
    df <- data.frame(id = rownames(mat), mat, check.names = FALSE)
    colnames(df)[1] <- idname
    write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
    message("Wrote ", path, " (", nrow(df), " rows x ", ncol(df) - 1, " samples)")
}

write_mat(round(txi$counts, 3),    snakemake@output[["gene_counts"]],  "gene_id")
write_mat(round(txi$abundance, 3), snakemake@output[["gene_tpm"]],     "gene_id")
write_mat(round(txi$length, 3),    snakemake@output[["gene_lengths"]], "gene_id")
write_mat(round(txi_tx$counts, 3), snakemake@output[["tx_counts"]],    "transcript_id")

zero_genes <- sum(rowSums(txi$counts) == 0)
message(sprintf("Genes with zero counts across all samples: %d / %d",
                zero_genes, nrow(txi$counts)))
if (zero_genes == nrow(txi$counts)) {
    stop("Every gene has zero counts. This usually means the reads did not map - ",
         "check the Salmon mapping rate in logs/salmon_quant/.")
}
