# Outputs

```
<outdir>/
├── 00_qc/
│   ├── multiqc_report.html        interactive QC across all samples and stages
│   ├── fastp/<sample>.fastp.{json,html}
│   └── fastqc/<sample>_R<n>.trimmed_fastqc.{html,zip}
├── 00_index/
│   ├── salmon/ or star/           built index (skipped if you supply a prebuilt one)
│   └── tx2gene.tsv                transcript → gene → gene_name, parsed from the GTF
├── 02_trimmed/                    adapter/quality-trimmed reads
├── 03_quant/
│   ├── salmon/<sample>/quant.sf   per-sample transcript quantification
│   ├── star/<sample>.*.bam        genomic + transcriptome BAM (star_salmon only)
│   ├── gene_counts.tsv            gene × sample counts  ← DESeq2 input
│   ├── gene_tpm.tsv               gene × sample TPM     ← for cross-gene comparison
│   ├── gene_lengths.tsv           effective lengths used in the correction
│   └── transcript_counts.tsv      transcript-level counts
├── 04_de/
│   ├── coldata.tsv                the sample metadata actually used in the model
│   ├── dds.rds                    DESeqDataSet - reload in R for custom analysis
│   ├── normalized_counts.tsv      median-of-ratios normalised counts
│   ├── vst.tsv                    variance-stabilised values (PCA/clustering/heatmaps)
│   ├── coefficients.txt           model coefficient names available to `results()`
│   └── <contrast>/
│       ├── results.tsv            all tested genes
│       └── significant.tsv        genes passing alpha (and lfc_threshold)
├── 05_report/
│   ├── rnaseq_report.html/.pdf    the client-facing report
│   └── summary_tables/*.tsv       the table behind every figure
└── logs/                          one log per rule
```

## Which file do I actually want?

- **Send to a collaborator/client** → `05_report/rnaseq_report.html` (self-contained;
  a single file, no dependencies) or the PDF.
- **Gene list for follow-up** → `04_de/<contrast>/significant.tsv`.
- **Never filter on `significant.tsv` alone when writing up** → `results.tsv` has every
  gene, including those that failed significance, which is what you need to say
  "gene X was not changed" honestly.
- **Plot expression of a specific gene** → `03_quant/gene_tpm.tsv` (comparable across
  genes) or `04_de/normalized_counts.tsv` (comparable across samples for one gene).
- **Custom modelling in R** → `04_de/dds.rds`:
  ```r
  dds <- readRDS("04_de/dds.rds")
  res <- DESeq2::results(dds, contrast = c("condition", "treated", "control"))
  ```

## results.tsv columns

| Column | Meaning |
|---|---|
| `gene_id` / `gene_name` | identifier and symbol (from the GTF) |
| `baseMean` | mean normalised count across all samples — how much evidence there is |
| `log2FoldChange` | shrunken effect size, numerator relative to denominator |
| `lfcSE` | standard error of that estimate |
| `pvalue` | raw p-value |
| `padj` | Benjamini–Hochberg adjusted p-value — **use this**, not `pvalue` |

`padj` is `NA` for genes filtered out by DESeq2's independent filtering (typically very
low counts); that means "not tested", not "not significant".

Fold changes are shrunken with `ashr`, which pulls in estimates for low-count genes.
This makes the ranking far more reliable but means the reported log2FC is deliberately
conservative relative to a raw ratio of means.
