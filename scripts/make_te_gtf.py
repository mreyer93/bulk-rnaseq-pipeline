"""Build a TEtranscripts-compatible transposable-element GTF from a UCSC RepeatMasker table.

Why this exists rather than downloading a prebuilt file: the TE GTFs that TEtranscripts
documents are distributed by hand (a lab file-share and a Dropbox link), and the
file-share path that most tutorials cite now returns 404. A pipeline that depends on that
link is a pipeline that stops working without warning. UCSC publishes the RepeatMasker
track for every assembly it hosts, under a stable URL, so generating the annotation is
both reproducible and version-pinnable.

The output follows the convention TEtranscripts expects:

    chr1  rmsk  exon  3000001  3000156  .  +  .  gene_id "L1Md_T"; transcript_id "L1Md_T_dup1"; family_id "L1"; class_id "LINE";

  gene_id       repeat name (the subfamily, e.g. L1Md_T) - counts aggregate to this
  transcript_id unique per locus, so each insertion is addressable
  family_id     repFamily (e.g. L1)
  class_id      repClass (e.g. LINE)

Run standalone:
    python scripts/make_te_gtf.py --genome mm39 --out mm39_rmsk_TE.gtf
    python scripts/make_te_gtf.py --rmsk rmsk.txt.gz --out TE.gtf
    python scripts/make_te_gtf.py --selftest
"""

# NB: no `from __future__ import annotations` here, for the same reason as tx2gene.py -
# Snakemake's `script:` directive prepends a preamble and a __future__ import must be
# the first statement in the file.

import argparse
import gzip
import io
import os
import sys
import urllib.request

UCSC_RMSK_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/{genome}/database/rmsk.txt.gz"

# UCSC rmsk column order. Most assemblies carry a leading 'bin' column; a few do not,
# so the parser keys off the column count rather than assuming.
_COLS_WITH_BIN = ["bin", "swScore", "milliDiv", "milliDel", "milliIns", "genoName",
                  "genoStart", "genoEnd", "genoLeft", "strand", "repName", "repClass",
                  "repFamily", "repStart", "repEnd", "repLeft", "id"]
_COLS_NO_BIN = _COLS_WITH_BIN[1:]

# Classes kept by default: the mobile elements people mean by "transposable element".
# Simple repeats, low-complexity regions, satellites and the small-RNA classes are
# excluded - they are repetitive but not TEs, and including them inflates the count
# matrix with features nobody will interpret.
DEFAULT_CLASSES = ("LINE", "SINE", "LTR", "DNA", "Retroposon", "RC")

EXCLUDED_BY_DEFAULT = ("Simple_repeat", "Low_complexity", "Satellite", "rRNA", "tRNA",
                       "snRNA", "srpRNA", "scRNA", "RNA", "Unknown", "Other")


def _open(path):
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def _clean_class(value):
    """Strip RepeatMasker's trailing '?' uncertainty marker: 'DNA?' -> 'DNA'."""
    return value[:-1] if value.endswith("?") else value


def parse_rmsk(handle, classes=DEFAULT_CLASSES, keep_ambiguous=True, min_length=0):
    """Yield (chrom, start1, end, strand, rep_name, rep_family, rep_class) per repeat.

    start1 is 1-based inclusive: UCSC stores genoStart 0-based half-open, GTF is 1-based.
    """
    wanted = set(classes) if classes else None
    for line in handle:
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= len(_COLS_WITH_BIN):
            cols = dict(zip(_COLS_WITH_BIN, parts))
        elif len(parts) >= len(_COLS_NO_BIN):
            cols = dict(zip(_COLS_NO_BIN, parts))
        else:
            continue

        raw_class = cols["repClass"]
        if not keep_ambiguous and raw_class.endswith("?"):
            continue
        rep_class = _clean_class(raw_class)
        if wanted is not None and rep_class not in wanted:
            continue

        try:
            start0 = int(cols["genoStart"])
            end = int(cols["genoEnd"])
        except ValueError:
            continue
        if end - start0 < min_length:
            continue

        strand = cols["strand"]
        # RepeatMasker writes 'C' for the complement strand in some dumps.
        if strand == "C":
            strand = "-"
        if strand not in ("+", "-"):
            strand = "."

        yield (cols["genoName"], start0 + 1, end, strand,
               cols["repName"], _clean_class(cols["repFamily"]), rep_class)


def write_gtf(records, out_handle):
    """Write TEtranscripts-format GTF. Returns (n_written, n_subfamilies, class_counts)."""
    dup = {}
    class_counts = {}
    n = 0
    for chrom, start, end, strand, name, family, klass in records:
        dup[name] = dup.get(name, 0) + 1
        out_handle.write(
            '{}\trmsk\texon\t{}\t{}\t.\t{}\t.\t'
            'gene_id "{}"; transcript_id "{}_dup{}"; family_id "{}"; class_id "{}";\n'
            .format(chrom, start, end, strand, name, name, dup[name], family, klass))
        class_counts[klass] = class_counts.get(klass, 0) + 1
        n += 1
    return n, len(dup), class_counts


def build(rmsk_path, out_path, classes=DEFAULT_CLASSES, keep_ambiguous=True, min_length=0):
    with _open(rmsk_path) as fh:
        opener = gzip.open(out_path, "wt") if out_path.endswith(".gz") else open(out_path, "w")
        with opener as out:
            return write_gtf(parse_rmsk(fh, classes, keep_ambiguous, min_length), out)


def download_rmsk(genome, dest):
    url = UCSC_RMSK_URL.format(genome=genome)
    sys.stderr.write("downloading {}\n".format(url))
    urllib.request.urlretrieve(url, dest)
    return dest


def _selftest():
    """Round-trip a handful of synthetic rmsk rows covering the cases that bite."""
    rows = [
        # bin  sw  div del ins  chrom  start   end     left strand name      class         family
        ["585", "1", "1", "0", "0", "chr1", "3000000", "3000156", "0", "+", "L1Md_T", "LINE", "L1"],
        ["585", "1", "1", "0", "0", "chr1", "3000200", "3000300", "0", "C", "L1Md_T", "LINE", "L1"],
        ["585", "1", "1", "0", "0", "chr1", "3000400", "3000500", "0", "-", "IAPEz", "LTR", "ERVK"],
        ["585", "1", "1", "0", "0", "chr1", "3000600", "3000700", "0", "+", "polyA", "Simple_repeat", "Simple_repeat"],
        ["585", "1", "1", "0", "0", "chr1", "3000800", "3000900", "0", "+", "MER5A", "DNA?", "hAT-Charlie"],
    ]
    lines = ["\t".join(r + ["1", "2", "3", "99"]) + "\n" for r in rows]

    out = io.StringIO()
    n, subfam, counts = write_gtf(parse_rmsk(iter(lines)), out)
    body = out.getvalue()
    checks = []

    checks.append(("simple repeats excluded", "polyA" not in body))
    checks.append(("LINE kept", 'gene_id "L1Md_T"' in body))
    checks.append(("LTR kept", 'class_id "LTR"' in body))
    checks.append(("ambiguous DNA? kept and cleaned", 'class_id "DNA"' in body and "DNA?" not in body))
    checks.append(("0-based start converted to 1-based", "\t3000001\t3000156\t" in body))
    checks.append(("C strand normalised to -", '\t3000201\t3000300\t.\t-\t' in body))
    checks.append(("loci uniquified per subfamily",
                   'transcript_id "L1Md_T_dup1"' in body and 'transcript_id "L1Md_T_dup2"' in body))
    checks.append(("counts", n == 4 and subfam == 3))

    ok = True
    for label, passed in checks:
        print("  {}  {}".format("PASS" if passed else "FAIL", label))
        ok = ok and passed

    # Negative control: the class filter must actually exclude when asked.
    out2 = io.StringIO()
    n2, _, _ = write_gtf(parse_rmsk(iter(lines), classes=("LTR",)), out2)
    only_ltr = n2 == 1 and 'gene_id "IAPEz"' in out2.getvalue()
    print("  {}  class filter restricts to requested classes".format("PASS" if only_ltr else "FAIL"))
    ok = ok and only_ltr

    # Negative control: keep_ambiguous=False must drop the DNA? record.
    out3 = io.StringIO()
    n3, _, _ = write_gtf(parse_rmsk(iter(lines), keep_ambiguous=False), out3)
    dropped = "MER5A" not in out3.getvalue()
    print("  {}  keep_ambiguous=False drops 'DNA?'".format("PASS" if dropped else "FAIL"))
    ok = ok and dropped

    print("\n{}".format("SELFTEST PASSED" if ok else "SELFTEST FAILED"))
    return 0 if ok else 1


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--genome", help="UCSC assembly to download rmsk for, e.g. mm39, hg38")
    p.add_argument("--rmsk", help="local UCSC rmsk table (.txt or .txt.gz)")
    p.add_argument("--out", help="output GTF (.gtf or .gtf.gz)")
    p.add_argument("--classes", default=",".join(DEFAULT_CLASSES),
                   help="comma-separated repClass values to keep, or 'all' "
                        "(default: %(default)s)")
    p.add_argument("--drop-ambiguous", action="store_true",
                   help="drop repeats whose class carries RepeatMasker's '?' marker")
    p.add_argument("--min-length", type=int, default=0,
                   help="skip repeats shorter than this many bp (default: 0)")
    p.add_argument("--selftest", action="store_true", help="run built-in checks and exit")
    args = p.parse_args(argv)

    if args.selftest:
        return _selftest()
    if not args.out:
        p.error("--out is required")
    if not args.rmsk and not args.genome:
        p.error("give either --rmsk or --genome")

    rmsk = args.rmsk
    tmp = None
    if not rmsk:
        tmp = args.out + ".rmsk.txt.gz"
        rmsk = download_rmsk(args.genome, tmp)

    classes = None if args.classes.strip().lower() == "all" else tuple(
        c.strip() for c in args.classes.split(",") if c.strip())

    n, subfam, counts = build(rmsk, args.out, classes, not args.drop_ambiguous,
                              args.min_length)
    sys.stderr.write("wrote {}: {:,} loci across {:,} subfamilies\n".format(args.out, n, subfam))
    for k in sorted(counts, key=lambda x: -counts[x]):
        sys.stderr.write("  {:<12s} {:>10,}\n".format(k, counts[k]))
    if n == 0:
        sys.stderr.write("ERROR: no repeats written - check --classes against the "
                         "repClass values in your rmsk table\n")
        return 1
    if tmp and os.path.exists(tmp):
        os.remove(tmp)
    return 0


if __name__ == "__main__":
    sys.exit(main())
