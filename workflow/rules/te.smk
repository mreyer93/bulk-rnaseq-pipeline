"""Transposable element and endogenous retrovirus quantification.

Off by default. Enable with `te_analysis: {enabled: true, ...}` in the config.

Why this needs its own alignment rather than reusing the star_salmon BAM: TE
quantification depends on multi-mapping reads, because young high-copy elements
(L1Md and IAP in mouse, L1HS and HERV-K in human) are near-identical across dozens to
thousands of loci. STAR's default `--outFilterMultimapNmax 10` discards a read mapping
to more than ten places, which is precisely the signal a TE analysis is looking for.
TEtranscripts' authors align at 100. So a BAM produced for gene-level work has already
thrown away what this analysis needs, and that includes essentially every BAM a
sequencing core will hand over.

Counting is TEcount (part of TEtranscripts), which assigns multi-mapping reads to TE
subfamilies by expectation-maximisation and emits genes and TEs in one table, so the
two can go into a single DESeq2 model and share one size factor.

Config:

    te_analysis:
      enabled: true
      te_gtf: "references/mm39_TE.gtf"    # from scripts/make_te_gtf.py
      multimap_max: 100                    # both STAR multimap limits
      mode: "multi"                        # TEcount --mode: multi | uniq
      exclude_bed: "references/blacklist.bed"   # optional, see below

`exclude_bed` exists because of a failure mode worth naming. Repeat-dense, gene-free
regions that are duplicated in the assembly act as sinks: reads pile up there carrying
low mapping quality, and a TE counter reports them as a large LINE or LTR signal that is
an artefact of the reference rather than biology. In one mouse study a single 10 kb
window took 7% of a library, every read at STAR MAPQ 3 (exactly two loci). Excluding
such regions is a judgement call, so it is explicit and logged rather than automatic.
"""

# TE_CFG and TE_ENABLED come from rules/common.smk, which parses them before
# rules/quantify.smk so that the STAR index rule knows it is needed.

if TE_ENABLED:
    TE_GTF = TE_CFG.get("te_gtf")
    if not TE_GTF:
        raise WorkflowError(
            "te_analysis.enabled is true but te_analysis.te_gtf is not set. Build one "
            "with:  python scripts/make_te_gtf.py --genome mm39 --out references/mm39_TE.gtf"
        )
    if not os.path.exists(TE_GTF):
        raise WorkflowError(f"te_analysis.te_gtf does not exist: {TE_GTF}")

    TE_MULTIMAP = int(TE_CFG.get("multimap_max", 100))
    TE_MODE = TE_CFG.get("mode", "multi")
    if TE_MODE not in ("multi", "uniq"):
        raise WorkflowError(f"te_analysis.mode must be 'multi' or 'uniq', got {TE_MODE!r}")
    TE_EXCLUDE = TE_CFG.get("exclude_bed", "")
    if TE_EXCLUDE and not os.path.exists(TE_EXCLUDE):
        raise WorkflowError(f"te_analysis.exclude_bed does not exist: {TE_EXCLUDE}")

    TE_DIR = join(OUTDIR, "06_te")

    rule te_star_align:
        """STAR with multi-mapping retained, for TE counting only.

        Deliberately separate from rule star_align: same index, different filters, and
        mixing the two would silently give whichever ran first.
        """
        input:
            reads = quant_input,
            index = STAR_INDEX,
        output:
            bam = join(TE_DIR, "bam", "{sample}.Aligned.sortedByCoord.out.bam"),
            log_final = join(TE_DIR, "bam", "{sample}.Log.final.out"),
        log: join(OUTDIR, "logs", "te_star_align", "{sample}.log")
        threads: config["threads"]["align"]
        params:
            prefix = join(TE_DIR, "bam", "{sample}."),
            nmax = TE_MULTIMAP,
            # TEtranscripts asks for winAnchorMultimapNmax to match outFilterMultimapNmax,
            # otherwise the anchor search caps the multimapping the filter would allow.
            extra = TE_CFG.get("star_extra", ""),
        conda: "../../envs/environment.yml"
        shell:
            """
            STAR --genomeDir {input.index} \
                --readFilesIn {input.reads} \
                --readFilesCommand zcat \
                --outFileNamePrefix {params.prefix} \
                --outSAMtype BAM SortedByCoordinate \
                --outFilterMultimapNmax {params.nmax} \
                --winAnchorMultimapNmax {params.nmax} \
                --outSAMattributes NH HI AS nM \
                --runThreadN {threads} {params.extra} > {log} 2>&1
            """

    rule te_filter_bam:
        """Drop alignments overlapping excluded regions, if a blacklist was given.

        A no-op passthrough when exclude_bed is unset, so the DAG shape does not change
        and the decision stays visible in one place.
        """
        input:
            bam = rules.te_star_align.output.bam,
        output:
            bam = join(TE_DIR, "bam_filtered", "{sample}.bam"),
        log: join(OUTDIR, "logs", "te_filter_bam", "{sample}.log")
        params:
            exclude = TE_EXCLUDE,
        conda: "../../envs/environment.yml"
        shell:
            """
            if [ -n "{params.exclude}" ]; then
                samtools view -b -L {params.exclude} -U {output.bam} \
                    -o /dev/null {input.bam} 2> {log}
            else
                ln -sf "$(cd "$(dirname {input.bam})" && pwd)/$(basename {input.bam})" \
                    {output.bam} 2> {log}
            fi
            """

    rule tecount:
        """Genes and TEs counted together, so both share one library size factor."""
        input:
            bam = rules.te_filter_bam.output.bam,
            gtf = config["reference"]["gtf"],
            te_gtf = TE_GTF,
        output:
            counts = join(TE_DIR, "counts", "{sample}.cntTable"),
        log: join(OUTDIR, "logs", "tecount", "{sample}.log")
        params:
            prefix = join(TE_DIR, "counts", "{sample}"),
            mode = TE_MODE,
            stranded = TE_CFG.get("stranded", "reverse"),
        conda: "../../envs/te.yml"
        shell:
            """
            TEcount --BAM {input.bam} \
                --GTF {input.gtf} \
                --TE {input.te_gtf} \
                --mode {params.mode} \
                --stranded {params.stranded} \
                --project {params.prefix} \
                --sortByPos > {log} 2>&1
            """

    rule te_matrix:
        """Aggregate per-sample cntTables into gene, TE and combined matrices."""
        input:
            counts = expand(join(TE_DIR, "counts", "{sample}.cntTable"),
                            sample=SAMPLE_NAMES),
            te_gtf = TE_GTF,
        output:
            combined = join(TE_DIR, "counts_combined.tsv"),
            te_only = join(TE_DIR, "te_counts.tsv"),
            gene_only = join(TE_DIR, "gene_counts.tsv"),
            summary = join(TE_DIR, "te_summary.tsv"),
        log: join(OUTDIR, "logs", "te_matrix.log")
        params:
            samples = SAMPLE_NAMES,
        conda: "../../envs/environment.yml"
        script: "../../scripts/te_matrix.py"
