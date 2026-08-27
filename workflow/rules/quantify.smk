"""Quantification. Two interchangeable paths, chosen by config['quantifier']:

  salmon       selective alignment straight from FASTQ against a transcriptome index.
               ~4-8 GB RAM for human. This is the local/laptop path.

  star_salmon  STAR genomic alignment, then Salmon quantification of the resulting
               transcriptome BAM. This is nf-core/rnaseq's default and gives you a
               genomic BAM for IGV/coverage work, but STAR needs ~30 GB RAM for a
               human index - hence the split. This is the cloud/server path.

Both converge on the same per-sample quant.sf, so everything downstream is identical.
"""


rule tx2gene:
    """Transcript -> gene mapping used by tximport, parsed from the GTF.

    Parsed here rather than assumed, because transcript IDs in a Salmon index must match
    the GTF exactly; deriving both from the same annotation avoids silent mismatches.
    """
    input: gtf = config["reference"].get("gtf", "")
    output: join(OUTDIR, "00_index", "tx2gene.tsv")
    log: join(OUTDIR, "logs", "tx2gene.log")
    conda: "../../envs/environment.yml"
    script: "../../scripts/tx2gene.py"


if QUANTIFIER == "salmon":

    if BUILD_SALMON_INDEX:
        rule salmon_index:
            input:
                txome = config["reference"].get("transcriptome_fasta", ""),
            output: directory(SALMON_INDEX)
            log: join(OUTDIR, "logs", "salmon_index.log")
            threads: config["threads"]["index"]
            params:
                kmer = config.get("salmon_kmer", 31),
                gencode = "--gencode" if config.get("gencode", False) else "",
            conda: "../../envs/environment.yml"
            shell:
                """
                salmon index \
                    -t {input.txome} \
                    -i {output} \
                    -k {params.kmer} \
                    -p {threads} {params.gencode} > {log} 2>&1
                """

    def _salmon_reads_arg(wildcards, input):
        if is_paired(wildcards.sample):
            return f"-1 {input.reads[0]} -2 {input.reads[1]}"
        return f"-r {input.reads[0]}"

    rule salmon_quant:
        input:
            reads = quant_input,
            index = SALMON_INDEX,
        output:
            quant = join(OUTDIR, "03_quant", "salmon", "{sample}", "quant.sf"),
        log: join(OUTDIR, "logs", "salmon_quant", "{sample}.log")
        threads: config["threads"]["quant"]
        params:
            reads_arg = _salmon_reads_arg,
            outdir = join(OUTDIR, "03_quant", "salmon", "{sample}"),
            libtype = config.get("salmon_libtype", "A"),  # A = auto-detect strandedness
            extra = config.get("salmon_extra", "--validateMappings --seqBias --gcBias"),
        conda: "../../envs/environment.yml"
        shell:
            """
            salmon quant \
                -i {input.index} \
                -l {params.libtype} \
                {params.reads_arg} \
                -o {params.outdir} \
                -p {threads} {params.extra} > {log} 2>&1
            """


elif QUANTIFIER == "star_salmon":

    if BUILD_STAR_INDEX:
        rule star_index:
            input:
                genome = config["reference"].get("genome_fasta", ""),
                gtf = config["reference"].get("gtf", ""),
            output: directory(STAR_INDEX)
            log: join(OUTDIR, "logs", "star_index.log")
            threads: config["threads"]["index"]
            params:
                sjdb = config.get("star_sjdb_overhang", 100),
                # Small genomes need genomeSAindexNbases reduced or STAR silently
                # produces a broken index; nf-core applies the same correction.
                extra = config.get("star_index_extra", ""),
                limit_ram = config.get("star_limit_ram", 0),
            conda: "../../envs/environment.yml"
            shell:
                """
                set -euo pipefail
                mkdir -p {output}
                GENOME_SIZE=$(grep -v '^>' {input.genome} | tr -d '\\n' | wc -c)
                # STAR's own recommendation: min(14, log2(GenomeLength)/2 - 1)
                NBASES=$(python3 -c "import math;print(min(14, int(math.log2($GENOME_SIZE)/2 - 1)))")
                echo "genome size ${{GENOME_SIZE}} -> --genomeSAindexNbases ${{NBASES}}" > {log}
                RAMARG=""
                if [ "{params.limit_ram}" != "0" ]; then
                    RAMARG="--limitGenomeGenerateRAM {params.limit_ram}"
                fi
                STAR --runMode genomeGenerate \
                    --genomeDir {output} \
                    --genomeFastaFiles {input.genome} \
                    --sjdbGTFfile {input.gtf} \
                    --sjdbOverhang {params.sjdb} \
                    --genomeSAindexNbases ${{NBASES}} \
                    --runThreadN {threads} \
                    ${{RAMARG}} {params.extra} >> {log} 2>&1
                """

    rule star_align:
        input:
            reads = quant_input,
            index = STAR_INDEX,
        output:
            bam = join(OUTDIR, "03_quant", "star", "{sample}.Aligned.toTranscriptome.out.bam"),
            genome_bam = join(OUTDIR, "03_quant", "star", "{sample}.Aligned.sortedByCoord.out.bam"),
            log_final = join(OUTDIR, "03_quant", "star", "{sample}.Log.final.out"),
        log: join(OUTDIR, "logs", "star_align", "{sample}.log")
        threads: config["threads"]["align"]
        params:
            prefix = join(OUTDIR, "03_quant", "star", "{sample}."),
            extra = config.get("star_align_extra", ""),
        conda: "../../envs/environment.yml"
        shell:
            """
            STAR --genomeDir {input.index} \
                --readFilesIn {input.reads} \
                --readFilesCommand zcat \
                --outFileNamePrefix {params.prefix} \
                --outSAMtype BAM SortedByCoordinate \
                --quantMode TranscriptomeSAM \
                --runThreadN {threads} {params.extra} > {log} 2>&1
            """

    rule salmon_quant_bam:
        """Salmon in alignment mode against the STAR transcriptome BAM."""
        input:
            bam = join(OUTDIR, "03_quant", "star", "{sample}.Aligned.toTranscriptome.out.bam"),
            txome = config["reference"].get("transcriptome_fasta", ""),
        output:
            quant = join(OUTDIR, "03_quant", "salmon", "{sample}", "quant.sf"),
        log: join(OUTDIR, "logs", "salmon_quant", "{sample}.log")
        threads: config["threads"]["quant"]
        params:
            outdir = join(OUTDIR, "03_quant", "salmon", "{sample}"),
            libtype = config.get("salmon_libtype", "A"),
        conda: "../../envs/environment.yml"
        shell:
            """
            salmon quant \
                -t {input.txome} \
                -l {params.libtype} \
                -a {input.bam} \
                -o {params.outdir} \
                -p {threads} > {log} 2>&1
            """


rule tximport:
    """Aggregate per-sample Salmon output to gene level (and keep transcript level).

    Uses tximport's lengthScaledTPM so the gene counts are appropriate as DESeq2 input
    while remaining comparable across samples.
    """
    input:
        quants = expand(join(OUTDIR, "03_quant", "salmon", "{sample}", "quant.sf"),
                        sample=SAMPLE_NAMES),
        tx2gene = join(OUTDIR, "00_index", "tx2gene.tsv"),
    output:
        gene_counts = join(OUTDIR, "03_quant", "gene_counts.tsv"),
        gene_tpm = join(OUTDIR, "03_quant", "gene_tpm.tsv"),
        gene_lengths = join(OUTDIR, "03_quant", "gene_lengths.tsv"),
        tx_counts = join(OUTDIR, "03_quant", "transcript_counts.tsv"),
    log: join(OUTDIR, "logs", "tximport.log")
    params:
        samples = SAMPLE_NAMES,
        quant_dir = join(OUTDIR, "03_quant", "salmon"),
    conda: "../../envs/r.yml"
    script: "../../scripts/tximport_counts.R"
