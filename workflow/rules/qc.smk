"""Read handling and QC: merge sequencing units, trim adapters, FastQC + MultiQC.

Single- and paired-end samples can coexist in one run. Which of the paired/single rules
applies to a given sample is decided by a wildcard constraint listing that sample set
explicitly - `ruleorder` cannot express this, because it resolves ambiguity globally
rather than per sample, and would force every sample down the paired path.
"""

import re as _re

PE_SAMPLES = [s for s in SAMPLE_NAMES if is_paired(s)]
SE_SAMPLES = [s for s in SAMPLE_NAMES if not is_paired(s)]


def _name_re(names):
    """Regex matching exactly this set of sample names (never-matching if empty)."""
    if not names:
        return r"(?!x)x"
    return "|".join(_re.escape(n) for n in names)


PE_RE = _name_re(PE_SAMPLES)
SE_RE = _name_re(SE_SAMPLES)


rule merge_fastq_r1:
    """Concatenate all sequencing units (lanes/runs) for a sample into one file.

    Concatenating gzip streams is valid gzip, so this works without recompressing.
    Samples with a single unit still pass through here, which keeps every downstream
    rule's input paths uniform.
    """
    input: lambda w: raw_fastqs(w, mate=1)
    output: temp(join(OUTDIR, "01_reads", "{sample}_R1.fastq.gz"))
    log: join(OUTDIR, "logs", "merge", "{sample}_R1.log")
    shell: "cat {input} > {output} 2> {log}"


rule merge_fastq_r2:
    wildcard_constraints: sample=PE_RE
    input: lambda w: raw_fastqs(w, mate=2)
    output: temp(join(OUTDIR, "01_reads", "{sample}_R2.fastq.gz"))
    log: join(OUTDIR, "logs", "merge", "{sample}_R2.log")
    shell: "cat {input} > {output} 2> {log}"


rule fastp_paired:
    """Adapter/quality trimming. fastp auto-detects adapters and emits a JSON that
    MultiQC picks up, so no adapter sequence needs to be configured."""
    wildcard_constraints: sample=PE_RE
    input:
        r1 = join(OUTDIR, "01_reads", "{sample}_R1.fastq.gz"),
        r2 = join(OUTDIR, "01_reads", "{sample}_R2.fastq.gz"),
    output:
        r1 = join(OUTDIR, "02_trimmed", "{sample}_R1.trimmed.fastq.gz"),
        r2 = join(OUTDIR, "02_trimmed", "{sample}_R2.trimmed.fastq.gz"),
        json = join(OUTDIR, "00_qc", "fastp", "{sample}.fastp.json"),
        html = join(OUTDIR, "00_qc", "fastp", "{sample}.fastp.html"),
    log: join(OUTDIR, "logs", "fastp", "{sample}.log")
    threads: config["threads"]["trim"]
    params:
        extra = config.get("fastp_extra", "")
    conda: "../../envs/environment.yml"
    shell:
        """
        fastp \
            -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            --json {output.json} --html {output.html} \
            --thread {threads} {params.extra} > {log} 2>&1
        """


rule fastp_single:
    wildcard_constraints: sample=SE_RE
    input:
        r1 = join(OUTDIR, "01_reads", "{sample}_R1.fastq.gz"),
    output:
        r1 = join(OUTDIR, "02_trimmed", "{sample}_R1.trimmed.fastq.gz"),
        json = join(OUTDIR, "00_qc", "fastp", "{sample}.fastp.json"),
        html = join(OUTDIR, "00_qc", "fastp", "{sample}.fastp.html"),
    log: join(OUTDIR, "logs", "fastp", "{sample}.log")
    threads: config["threads"]["trim"]
    params:
        extra = config.get("fastp_extra", "")
    conda: "../../envs/environment.yml"
    shell:
        """
        fastp \
            -i {input.r1} -o {output.r1} \
            --json {output.json} --html {output.html} \
            --thread {threads} {params.extra} > {log} 2>&1
        """


rule fastqc:
    input: join(OUTDIR, "02_trimmed", "{sample}_R{mate}.trimmed.fastq.gz")
    output:
        html = join(OUTDIR, "00_qc", "fastqc", "{sample}_R{mate}.trimmed_fastqc.html"),
        zip  = join(OUTDIR, "00_qc", "fastqc", "{sample}_R{mate}.trimmed_fastqc.zip"),
    log: join(OUTDIR, "logs", "fastqc", "{sample}_R{mate}.log")
    params: outdir = join(OUTDIR, "00_qc", "fastqc")
    conda: "../../envs/environment.yml"
    shell: "fastqc -o {params.outdir} -q {input} > {log} 2>&1"


def multiqc_inputs(wildcards):
    """Everything MultiQC should summarise, depending on how this run is configured."""
    files = []
    for s in SAMPLE_NAMES:
        files.append(join(OUTDIR, "00_qc", "fastp", f"{s}.fastp.json"))
        mates = [1, 2] if is_paired(s) else [1]
        for m in mates:
            files.append(join(OUTDIR, "00_qc", "fastqc",
                              f"{s}_R{m}.trimmed_fastqc.zip"))
        if QUANTIFIER == "salmon":
            files.append(join(OUTDIR, "03_quant", "salmon", s, "quant.sf"))
        else:
            files.append(join(OUTDIR, "03_quant", "star", f"{s}.Log.final.out"))
            files.append(join(OUTDIR, "03_quant", "salmon", s, "quant.sf"))
    return files


rule multiqc:
    """Single interactive QC page across trimming, read quality and quantification."""
    input: multiqc_inputs
    output: join(OUTDIR, "00_qc", "multiqc_report.html")
    log: join(OUTDIR, "logs", "multiqc.log")
    params:
        indir = OUTDIR,
        outdir = join(OUTDIR, "00_qc"),
    conda: "../../envs/environment.yml"
    shell:
        """
        multiqc --force --no-ansi \
            --outdir {params.outdir} \
            --filename multiqc_report.html \
            {params.indir} > {log} 2>&1
        """
