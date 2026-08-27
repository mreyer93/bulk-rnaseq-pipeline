"""Build a transcript -> gene mapping table from a GTF/GFF annotation.

Written to be tolerant of the annotation flavours people actually hand you: Ensembl,
GENCODE (whose transcript IDs carry version suffixes and whose FASTA headers are
pipe-delimited), RefSeq, and the minimal GTFs that come out of StringTie or a genome
browser. Emits columns: transcript_id, gene_id, gene_name.

Run standalone for a quick check:
    python scripts/tx2gene.py annotation.gtf out.tsv
"""

# NB: no `from __future__ import annotations` here. Snakemake's `script:` directive
# prepends its own preamble to this file before executing it, which would push a
# __future__ import off line 1 and raise SyntaxError. Keep type hints 3.8-compatible.

import gzip
import io
import os
import re
import sys

_ATTR_GTF = re.compile(r'(\S+)\s+"([^"]*)"')      # GTF:  key "value";
_ATTR_GFF = re.compile(r'([^=;]+)=([^;]*)')        # GFF3: key=value;


def _open(path):
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def _attrs(field: str) -> dict:
    d = dict(_ATTR_GTF.findall(field))
    if not d:
        d = {k.strip(): v.strip() for k, v in _ATTR_GFF.findall(field)}
    return d


def parse(gtf_path: str, strip_version: bool = False):
    """Return (rows, stats). rows = list of (transcript_id, gene_id, gene_name)."""
    seen = {}
    n_lines = n_feat = 0
    with _open(gtf_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            n_lines += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feature, attr_field = parts[2], parts[8]
            # Prefer transcript features, but fall back to exons: some minimal GTFs
            # (and a few RefSeq conversions) have no standalone transcript lines.
            if feature not in ("transcript", "exon", "mRNA"):
                continue
            a = _attrs(attr_field)
            tx = a.get("transcript_id") or a.get("transcript") or a.get("ID")
            if not tx:
                continue
            gene = (a.get("gene_id") or a.get("gene") or a.get("Parent") or tx)
            name = (a.get("gene_name") or a.get("gene_symbol") or a.get("Name") or gene)
            if strip_version:
                tx = tx.split(".")[0]
                gene = gene.split(".")[0]
            if tx not in seen:
                seen[tx] = (tx, gene, name)
                n_feat += 1
    return list(seen.values()), {"lines": n_lines, "transcripts": n_feat}


def write(rows, out_path):
    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("transcript_id\tgene_id\tgene_name\n")
        for tx, gene, name in rows:
            fh.write(f"{tx}\t{gene}\t{name}\n")


def run(gtf_path, out_path, strip_version=False, log=None):
    def _log(msg):
        print(msg)
        if log:
            with open(log, "a", encoding="utf-8") as fh:
                fh.write(msg + "\n")

    if not gtf_path or not os.path.exists(gtf_path):
        raise SystemExit(f"ERROR: annotation not found: {gtf_path!r}")

    rows, stats = parse(gtf_path, strip_version=strip_version)
    if not rows:
        raise SystemExit(
            f"ERROR: no transcript/gene pairs parsed from {gtf_path}.\n"
            f"Read {stats['lines']} feature lines but found no transcript_id attribute. "
            f"Is this really a GTF/GFF annotation?"
        )
    write(rows, out_path)
    n_genes = len({r[1] for r in rows})
    _log(f"Parsed {stats['transcripts']} transcripts across {n_genes} genes "
         f"from {gtf_path} -> {out_path}")
    return rows


# --- Snakemake entry point --------------------------------------------------------
if "snakemake" in globals():
    run(
        snakemake.input.gtf,                                    # noqa: F821
        snakemake.output[0],                                    # noqa: F821
        strip_version=snakemake.config.get("strip_tx_version", False),  # noqa: F821
        log=snakemake.log[0] if snakemake.log else None,        # noqa: F821
    )
elif __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__)
        raise SystemExit(2)
    run(sys.argv[1], sys.argv[2], strip_version="--strip-version" in sys.argv)
