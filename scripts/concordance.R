#!/usr/bin/env Rscript
# Cross-tissue concordance: genes moving the same way in two independent sets of contrasts.
#
# The motivating question is a peripheral-biomarker question. Something is measured in a
# tissue that cannot be sampled in a living subject (brain, tumour, bone marrow) and the
# goal is a gene whose change there is mirrored in a tissue that can be sampled (blood).
# That is not a standard differential expression output, and the obvious approach to it is
# wrong in a way that matters.
#
# The obvious approach is to threshold each tissue at FDR < 0.05 and intersect the two
# lists. Two problems. First, the intersection of two 5% FDR lists does not itself have a
# 5% false discovery rate; the error rate of the intersection is not controlled at any
# stated level, and it is usually reported as though it were. Second, thresholding throws
# away every gene that lands just outside the cut in one tissue while being strongly and
# consistently moved in both, which is exactly the population a biomarker is likely to sit
# in.
#
# Three complementary analyses are produced instead:
#
#   1. Intersection-union test. The hypothesis "changed in tissue A AND changed in
#      tissue B" is a conjunction, and the valid test statistic for a conjunction is the
#      MAXIMUM of the component p-values, not their combination. Fisher's or Stouffer's
#      method answers a different question ("changed in at least one"), and using either
#      here inflates the apparent evidence. max(p) is conservative and correct, and FDR is
#      then applied to it. This is the primary result.
#
#   2. Threshold-free effect-size agreement. Spearman and Pearson correlation of log2 fold
#      changes across all tested genes, which says whether the two tissues respond
#      similarly at all, independent of any cutoff.
#
#   3. Directional concordance. Of genes passing a lenient screen, what fraction move the
#      same way, against the 50% expected by chance.
#
# Usage, standalone:
#   Rscript scripts/concordance.R --a A1.tsv,A2.tsv --b B1.tsv,B2.tsv --outdir out/
#   Rscript scripts/concordance.R --selftest
#
# Input tables are DESeq2 results as written by scripts/deseq2.R: a gene id column plus
# log2FoldChange, pvalue and padj. Multiple files per side are combined by taking, per
# gene, the least significant p-value across that side's contrasts. That is deliberate:
# requiring a gene to hold up across every contrast within a tissue group is the stricter
# and more useful claim for a biomarker.

suppressPackageStartupMessages({
  library(stats)
})

# ---------------------------------------------------------------- helpers ---

read_de <- function(path) {
  df <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  id_col <- intersect(c("gene", "gene_id", "id", "feature", "X"), names(df))
  if (length(id_col) == 0) id_col <- names(df)[1] else id_col <- id_col[1]
  need <- c("log2FoldChange", "pvalue")
  missing <- setdiff(need, names(df))
  if (length(missing)) {
    stop(sprintf("%s is missing column(s): %s", path, paste(missing, collapse = ", ")))
  }
  data.frame(gene = as.character(df[[id_col]]),
             lfc = as.numeric(df$log2FoldChange),
             p = as.numeric(df$pvalue),
             stringsAsFactors = FALSE)
}

# Collapse several contrasts from one tissue group into a single per-gene summary.
# worst-case p (max) and mean lfc: a gene must hold up in every contrast in the group.
collapse_side <- function(paths) {
  tabs <- lapply(paths, read_de)
  genes <- Reduce(intersect, lapply(tabs, function(t) t$gene))
  if (length(genes) == 0) stop("no genes shared across the contrasts on one side")
  lfc <- sapply(tabs, function(t) t$lfc[match(genes, t$gene)])
  p   <- sapply(tabs, function(t) t$p[match(genes, t$gene)])
  if (is.null(dim(lfc))) { lfc <- matrix(lfc, ncol = 1); p <- matrix(p, ncol = 1) }
  # NA p-values are DESeq2's independent-filtering outcome; treat as no evidence.
  p[is.na(p)] <- 1
  data.frame(gene = genes,
             lfc = rowMeans(lfc, na.rm = TRUE),
             p = apply(p, 1, max),
             stringsAsFactors = FALSE)
}

concordance <- function(a, b, direction = "up", alpha = 0.05, screen_p = 0.05) {
  genes <- intersect(a$gene, b$gene)
  if (length(genes) == 0) stop("the two sides share no genes")
  A <- a[match(genes, a$gene), ]
  B <- b[match(genes, b$gene), ]

  # 1. Intersection-union test: the conjunction p-value is the maximum, not a combination.
  p_conj <- pmax(A$p, B$p)
  same_sign <- sign(A$lfc) == sign(B$lfc)
  dir_ok <- switch(direction,
                   up   = A$lfc > 0 & B$lfc > 0,
                   down = A$lfc < 0 & B$lfc < 0,
                   any  = same_sign)
  # A gene moving in opposite directions cannot satisfy a conjunction in one direction.
  p_conj[!dir_ok] <- 1
  padj_conj <- p.adjust(p_conj, method = "BH")

  res <- data.frame(gene = genes,
                    lfc_a = A$lfc, p_a = A$p,
                    lfc_b = B$lfc, p_b = B$p,
                    lfc_mean = (A$lfc + B$lfc) / 2,
                    same_direction = same_sign,
                    p_conjunction = p_conj,
                    padj_conjunction = padj_conj,
                    stringsAsFactors = FALSE)
  res <- res[order(res$padj_conjunction, -abs(res$lfc_mean)), ]

  # 2. Threshold-free agreement across all shared genes.
  ok <- is.finite(A$lfc) & is.finite(B$lfc)
  cors <- list(
    spearman = suppressWarnings(cor(A$lfc[ok], B$lfc[ok], method = "spearman")),
    pearson  = suppressWarnings(cor(A$lfc[ok], B$lfc[ok], method = "pearson")),
    n = sum(ok))

  # 3. Directional concordance among genes with some evidence on both sides, tested
  #    against the 50% expected if the two tissues were unrelated.
  screen <- A$p < screen_p & B$p < screen_p & ok
  n_screen <- sum(screen)
  n_same <- sum(same_sign[screen])
  bt <- if (n_screen > 0) binom.test(n_same, n_screen, p = 0.5) else NULL

  list(results = res, cors = cors,
       n_screen = n_screen, n_same = n_same,
       binom_p = if (is.null(bt)) NA_real_ else bt$p.value,
       n_hits = sum(res$padj_conjunction < alpha, na.rm = TRUE),
       alpha = alpha, direction = direction)
}

# ---------------------------------------------------------------- plotting ---

plot_concordance <- function(cc, path, label_a = "tissue A", label_b = "tissue B",
                             top_n = 12) {
  r <- cc$results
  ok <- is.finite(r$lfc_a) & is.finite(r$lfc_b)
  r <- r[ok, ]
  hit <- r$padj_conjunction < cc$alpha

  png(path, width = 1500, height = 1400, res = 200)
  on.exit(dev.off(), add = TRUE)
  par(mar = c(4.6, 4.6, 4.2, 1.2))
  lim <- quantile(abs(c(r$lfc_a, r$lfc_b)), 0.999, na.rm = TRUE)
  plot(r$lfc_a, r$lfc_b, pch = 16, cex = 0.35,
       col = adjustcolor("#9AA5AD", 0.45),
       xlim = c(-lim, lim), ylim = c(-lim, lim),
       xlab = sprintf("log2 fold change, %s", label_a),
       ylab = sprintf("log2 fold change, %s", label_b),
       main = "")
  abline(h = 0, v = 0, col = "#CCCCCC")
  abline(a = 0, b = 1, lty = 2, col = "#CCCCCC")
  if (any(hit)) {
    points(r$lfc_a[hit], r$lfc_b[hit], pch = 16, cex = 0.8, col = "#C1442E")
    lab <- head(which(hit), top_n)
    text(r$lfc_a[lab], r$lfc_b[lab], labels = r$gene[lab],
         pos = 4, cex = 0.55, col = "#7A2A1C", xpd = NA)
  }
  title(main = sprintf("Cross-tissue concordance: %d gene%s at FDR < %.2f",
                       cc$n_hits, ifelse(cc$n_hits == 1, "", "s"), cc$alpha),
        adj = 0, cex.main = 1.05, line = 2.6)
  title(main = sprintf("Spearman rho = %.3f over %s genes; %.0f%% of co-moving genes agree in direction",
                       cc$cors$spearman, format(cc$cors$n, big.mark = ","),
                       100 * cc$n_same / max(cc$n_screen, 1)),
        adj = 0, cex.main = 0.78, font.main = 1, line = 1.4, col.main = "#555555")
  invisible(path)
}

# ---------------------------------------------------------------- selftest ---

.selftest <- function() {
  set.seed(1)
  n <- 3000
  genes <- sprintf("G%04d", seq_len(n))
  # 100 genes genuinely co-regulated up in both; the rest noise.
  lfc_a <- rnorm(n, 0, 0.3); lfc_b <- rnorm(n, 0, 0.3)
  p_a <- runif(n); p_b <- runif(n)
  true_idx <- 1:100
  lfc_a[true_idx] <- runif(100, 1.2, 3); lfc_b[true_idx] <- runif(100, 1.0, 2.5)
  p_a[true_idx] <- runif(100, 0, 1e-6);  p_b[true_idx] <- runif(100, 0, 1e-6)
  # 60 genes strongly moved in OPPOSITE directions: must never be called concordant-up.
  opp_idx <- 101:160
  lfc_a[opp_idx] <- runif(60, 1.5, 3); lfc_b[opp_idx] <- -runif(60, 1.5, 3)
  p_a[opp_idx] <- runif(60, 0, 1e-8);  p_b[opp_idx] <- runif(60, 0, 1e-8)

  A <- data.frame(gene = genes, lfc = lfc_a, p = p_a, stringsAsFactors = FALSE)
  B <- data.frame(gene = genes, lfc = lfc_b, p = p_b, stringsAsFactors = FALSE)
  cc <- concordance(A, B, direction = "up", alpha = 0.05)

  hits <- cc$results$gene[cc$results$padj_conjunction < 0.05]
  checks <- list(
    c("recovers co-regulated genes", sum(hits %in% genes[true_idx]) >= 95),
    c("excludes opposite-direction genes", !any(hits %in% genes[opp_idx])),
    # The max-p invariant holds for every gene that passes the direction filter.
    # Genes failing it are set to p = 1 by construction, which is the point of the
    # filter, so they are excluded from this check rather than treated as failures.
    c("conjunction p is the max of the two, where direction allows", {
      r <- cc$results
      keep <- r$lfc_a > 0 & r$lfc_b > 0
      isTRUE(all.equal(r$p_conjunction[keep], pmax(r$p_a[keep], r$p_b[keep])))
    }),
    c("direction-failing genes are neutralised to p = 1", {
      r <- cc$results
      drop <- !(r$lfc_a > 0 & r$lfc_b > 0)
      all(r$p_conjunction[drop] == 1)
    }),
    c("false positives controlled", sum(!(hits %in% genes[true_idx])) <= 5),
    c("direction flag honoured", all(cc$results$lfc_a[cc$results$padj_conjunction < 0.05] > 0))
  )

  # Negative control: an unrelated pair must yield essentially nothing.
  B2 <- data.frame(gene = genes, lfc = rnorm(n, 0, 0.3), p = runif(n),
                   stringsAsFactors = FALSE)
  cc2 <- concordance(A, B2, direction = "up", alpha = 0.05)
  checks[[length(checks) + 1]] <- c("null pair yields few hits", cc2$n_hits <= 5)

  # Negative control: Fisher's method would wrongly call the opposite-direction genes.
  # Demonstrates why the maximum is used rather than a combination.
  fisher_p <- pchisq(-2 * (log(p_a[opp_idx]) + log(p_b[opp_idx])), df = 4, lower.tail = FALSE)
  checks[[length(checks) + 1]] <- c("Fisher would have called the discordant genes",
                                    all(fisher_p < 0.05))

  ok <- TRUE
  for (ch in checks) {
    passed <- isTRUE(as.logical(ch[[2]]))
    cat(sprintf("  %s  %s\n", ifelse(passed, "PASS", "FAIL"), ch[[1]]))
    ok <- ok && passed
  }
  cat(sprintf("\n%s\n", ifelse(ok, "SELFTEST PASSED", "SELFTEST FAILED")))
  if (!ok) quit(status = 1)
  invisible(TRUE)
}

# ---------------------------------------------------------------- main ---

# Run the command-line interface only when this file is the script being executed.
# `interactive()` is not the right test: it is FALSE under Rscript whether this file is
# the entry point or was source()d by another script, so using it made the CLI fire (and
# fail on missing --a/--b) whenever another analysis tried to reuse these functions as a
# library. Comparing against --file= distinguishes the two cases.
.called_directly <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  length(f) == 1 && basename(f) == "concordance.R"
})

if (.called_directly) {
  args <- commandArgs(trailingOnly = TRUE)
  if ("--selftest" %in% args) {
    .selftest()
  } else if (exists("snakemake")) {
    a <- collapse_side(as.character(snakemake@params[["group_a"]]))
    b <- collapse_side(as.character(snakemake@params[["group_b"]]))
    cc <- concordance(a, b,
                      direction = snakemake@params[["direction"]],
                      alpha = as.numeric(snakemake@params[["alpha"]]))
    write.table(cc$results, snakemake@output[["table"]], sep = "\t",
                quote = FALSE, row.names = FALSE)
    plot_concordance(cc, snakemake@output[["figure"]],
                     snakemake@params[["label_a"]], snakemake@params[["label_b"]])
  } else {
    get_arg <- function(flag, default = NULL) {
      i <- match(flag, args)
      if (is.na(i) || i == length(args)) default else args[i + 1]
    }
    a_files <- strsplit(get_arg("--a", ""), ",")[[1]]
    b_files <- strsplit(get_arg("--b", ""), ",")[[1]]
    outdir <- get_arg("--outdir", ".")
    if (!length(a_files) || !length(b_files)) {
      stop("give --a and --b as comma-separated DESeq2 results tables")
    }
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    cc <- concordance(collapse_side(a_files), collapse_side(b_files),
                      direction = get_arg("--direction", "up"),
                      alpha = as.numeric(get_arg("--alpha", "0.05")))
    write.table(cc$results, file.path(outdir, "concordance.tsv"), sep = "\t",
                quote = FALSE, row.names = FALSE)
    plot_concordance(cc, file.path(outdir, "concordance.png"),
                     get_arg("--label-a", "tissue A"), get_arg("--label-b", "tissue B"))
    cat(sprintf("concordance: %d genes at FDR < %.2f; Spearman rho %.3f over %s genes\n",
                cc$n_hits, cc$alpha, cc$cors$spearman, format(cc$cors$n, big.mark = ",")))
  }
}
