#!/usr/bin/env Rscript
# Regenerates the figures shown in example/README.md from a completed pipeline run.
# Run after ./test/run_test.sh:   Rscript example/make_figures.R test/data/results
suppressPackageStartupMessages({ library(ggplot2) })

res <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(res)) res <- "test/data/results"
figdir <- "example/figures"; outdir <- "example/results"
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

PAL <- c("#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#937860")
th <- theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey30", size = 10))
sv <- function(p, name, w = 8, h = 5) {
  ggsave(file.path(figdir, paste0(name, ".png")), p, width = w, height = h, dpi = 150, bg = "white")
  cat("  wrote", name, "\n")
}

vst <- read.delim(file.path(res, "04_de/vst.tsv"), row.names = 1, check.names = FALSE)
cold <- read.delim(file.path(res, "04_de/coldata.tsv"), row.names = 1, check.names = FALSE)
cts <- read.delim(file.path(res, "03_quant/gene_counts.tsv"), row.names = 1, check.names = FALSE)

# --- 1. PCA -----------------------------------------------------------------
m <- as.matrix(vst); v <- apply(m, 1, var)
keep <- head(order(v, decreasing = TRUE), min(500, sum(v > 0)))
pca <- prcomp(t(m[keep, ]), scale. = FALSE)
pv <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pd <- data.frame(sample = rownames(pca$x), PC1 = pca$x[,1], PC2 = pca$x[,2],
                 condition = cold[rownames(pca$x), "condition"])
# replicate number is enough on the label; the condition is already the colour
pd$lab <- sub("^.*_(REP[0-9]+)$", "\\1", pd$sample)
xr <- range(pd$PC1); yr <- range(pd$PC2)
sv(ggplot(pd, aes(PC1, PC2, colour = condition)) +
     geom_point(size = 4) +
     geom_text(aes(label = lab), vjust = -1.3, size = 3.4, show.legend = FALSE) +
     scale_colour_manual(values = PAL) +
     expand_limits(x = xr + c(-1, 1) * diff(xr) * 0.12,
                   y = yr + c(-1, 1) * diff(yr) * 0.22) +
     labs(title = "Samples separate by condition",
          subtitle = "PCA on the 500 most variable genes, variance-stabilised",
          x = sprintf("PC1 (%.1f%%)", pv[1]), y = sprintf("PC2 (%.1f%%)", pv[2])) + th,
   "pca")

# --- 2. Library size / mapping ------------------------------------------------
ls_df <- data.frame(sample = colnames(cts), counts = colSums(cts),
                    genes = colSums(cts > 0))
ls_df$condition <- cold[ls_df$sample, "condition"]
write.table(ls_df, file.path(outdir, "library_sizes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
sv(ggplot(ls_df, aes(reorder(sample, counts), counts, fill = condition)) +
     geom_col() + coord_flip() + scale_fill_manual(values = PAL) +
     labs(title = "Library size by sample", subtitle = "assigned counts after quantification",
          x = NULL, y = "Counts") + th, "library_sizes", 8, 4)

# --- 3. Volcano + 4. heatmap, per contrast -----------------------------------
contrasts <- list.dirs(file.path(res, "04_de"), recursive = FALSE, full.names = FALSE)
contrasts <- contrasts[grepl("_vs_", contrasts)]
for (cn in contrasts) {
  d <- read.delim(file.path(res, "04_de", cn, "results.tsv"), check.names = FALSE)
  d <- d[!is.na(d$padj), ]
  if (!nrow(d)) next
  d$sig <- d$padj < 0.05
  d$dir <- ifelse(!d$sig, "n.s.", ifelse(d$log2FoldChange > 0, "up", "down"))
  d$nlp <- -log10(pmax(d$padj, .Machine$double.xmin))
  lab <- head(d[order(d$padj), ], 8)
  sv(ggplot(d, aes(log2FoldChange, nlp, colour = dir)) +
       geom_point(alpha = 0.75, size = 2) +
       geom_text(data = lab, aes(label = gene_name), vjust = -0.8, size = 2.9, show.legend = FALSE) +
       scale_colour_manual(values = c(up = PAL[4], down = PAL[1], `n.s.` = "grey75"), name = NULL) +
       geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
       labs(title = paste("Differential expression:", gsub("_", " ", cn)),
            subtitle = sprintf("%d of %d genes significant at adjusted p < 0.05",
                               sum(d$sig), nrow(d)),
            x = "log2 fold change", y = "-log10 adjusted p") + th,
     paste0("volcano_", cn))
  write.table(head(d[order(d$padj), c("gene_id","gene_name","baseMean","log2FoldChange","padj")], 15),
              file.path(outdir, paste0("top_genes_", cn, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
}

# --- 5. Heatmap of top genes (first contrast) --------------------------------
cn <- contrasts[1]
d <- read.delim(file.path(res, "04_de", cn, "results.tsv"), check.names = FALSE)
d <- d[!is.na(d$padj), ]; d <- d[order(d$padj), ]
g <- head(intersect(d$gene_id, rownames(vst)), 25)
mm <- as.matrix(vst[g, ]); mm <- t(scale(t(mm), center = TRUE, scale = FALSE))
lbl <- d$gene_name[match(g, d$gene_id)]; lbl[is.na(lbl) | lbl == ""] <- g[is.na(lbl) | lbl == ""]
hd <- expand.grid(gene = lbl, sample = colnames(mm), stringsAsFactors = FALSE)
hd$value <- as.vector(mm); hd$gene <- factor(hd$gene, levels = rev(lbl))
sv(ggplot(hd, aes(sample, gene, fill = value)) + geom_tile() +
     scale_fill_gradient2(low = PAL[1], mid = "white", high = PAL[4], midpoint = 0,
                          name = "centred\nexpression") +
     labs(title = "Top differentially expressed genes",
          subtitle = paste("ranked by adjusted p,", gsub("_", " ", cn)), x = NULL, y = NULL) +
     th + theme(axis.text.x = element_text(angle = 40, hjust = 1)), "heatmap_top_genes", 8, 6.5)

cat("done\n")
