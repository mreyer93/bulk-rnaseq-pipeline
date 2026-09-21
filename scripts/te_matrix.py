"""Aggregate per-sample TEcount tables into gene, TE and combined count matrices.

TEcount emits one `.cntTable` per sample with a single count column. Gene rows carry the
plain gene identifier; TE rows carry a compound `name:family:class` identifier, for
example `L1Md_T:L1:LINE`. That compound form is how the two are told apart here, and it
is also what lets the TE rows be rolled up by family and class without a second
annotation lookup.

Genes and TEs are kept in one combined matrix as well as split, because the right way to
normalise them is together: a single DESeq2 size factor across both means a global shift
in TE expression is measured against the same library scaling as the genes, rather than
against a separately normalised TE-only library that would partly absorb the very effect
being looked for.

Run standalone:
    python scripts/te_matrix.py --selftest
"""

# NB: no `from __future__ import annotations` - Snakemake's script: directive prepends a
# preamble, which would stop a __future__ import being the first statement.

import os
import sys


def parse_cnt_table(path):
    """Return {feature_id: count} from one TEcount .cntTable."""
    counts = {}
    with open(path) as fh:
        header = fh.readline()          # 'gene/TE<TAB><sample>'
        if header and "\t" not in header:
            raise ValueError("%s does not look like a TEcount table" % path)
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            try:
                counts[parts[0]] = int(float(parts[1]))
            except ValueError:
                continue
    return counts


def is_te(feature_id):
    """TE rows are 'name:family:class'; gene rows are bare identifiers.

    Checked on the colon count rather than a suffix match so that gene IDs containing a
    colon (rare, but they exist in some annotations) are not silently reclassified.
    """
    return feature_id.count(":") == 2


def split_te_id(feature_id):
    name, family, klass = feature_id.split(":")
    return name, family, klass


def build(sample_files, samples):
    """sample_files aligned with samples. Returns (features, matrix) with matrix[f][s]."""
    per_sample = {}
    for s, path in zip(samples, sample_files):
        per_sample[s] = parse_cnt_table(path)

    features = []
    seen = set()
    for s in samples:
        for f in per_sample[s]:
            if f not in seen:
                seen.add(f)
                features.append(f)

    matrix = {f: {s: per_sample[s].get(f, 0) for s in samples} for f in features}
    return features, matrix


def write_matrix(features, matrix, samples, path, header="feature"):
    with open(path, "w") as out:
        out.write(header + "\t" + "\t".join(samples) + "\n")
        for f in features:
            out.write(f + "\t" + "\t".join(str(matrix[f][s]) for s in samples) + "\n")
    return len(features)


def summarise_te(features, matrix, samples, path):
    """Per-sample TE totals by class and by family, plus the TE share of all counts."""
    by_class = {}
    by_family = {}
    te_total = dict((s, 0) for s in samples)
    gene_total = dict((s, 0) for s in samples)

    for f in features:
        if is_te(f):
            _, family, klass = split_te_id(f)
            for s in samples:
                v = matrix[f][s]
                by_class.setdefault(klass, dict((x, 0) for x in samples))[s] += v
                by_family.setdefault(family, dict((x, 0) for x in samples))[s] += v
                te_total[s] += v
        else:
            for s in samples:
                gene_total[s] += matrix[f][s]

    with open(path, "w") as out:
        out.write("level\tcategory\t" + "\t".join(samples) + "\n")
        out.write("total\tgene\t" + "\t".join(str(gene_total[s]) for s in samples) + "\n")
        out.write("total\tTE\t" + "\t".join(str(te_total[s]) for s in samples) + "\n")
        out.write("fraction\tTE_share\t" + "\t".join(
            "%.5f" % (te_total[s] / max(te_total[s] + gene_total[s], 1)) for s in samples) + "\n")
        for klass in sorted(by_class, key=lambda k: -sum(by_class[k].values())):
            out.write("class\t%s\t" % klass +
                      "\t".join(str(by_class[klass][s]) for s in samples) + "\n")
        for family in sorted(by_family, key=lambda k: -sum(by_family[k].values()))[:60]:
            out.write("family\t%s\t" % family +
                      "\t".join(str(by_family[family][s]) for s in samples) + "\n")
    return te_total, gene_total


def run(sample_files, samples, out_combined, out_te, out_gene, out_summary):
    features, matrix = build(sample_files, samples)
    te_features = [f for f in features if is_te(f)]
    gene_features = [f for f in features if not is_te(f)]

    write_matrix(features, matrix, samples, out_combined)
    write_matrix(te_features, matrix, samples, out_te, header="TE")
    write_matrix(gene_features, matrix, samples, out_gene, header="gene")
    te_total, gene_total = summarise_te(features, matrix, samples, out_summary)

    sys.stderr.write(
        "te_matrix: %d samples, %d genes, %d TE subfamilies\n"
        % (len(samples), len(gene_features), len(te_features)))
    for s in samples:
        tot = te_total[s] + gene_total[s]
        sys.stderr.write("  %-24s TE share %.2f%%  (%s TE / %s total counts)\n"
                         % (s, 100.0 * te_total[s] / max(tot, 1),
                            format(te_total[s], ","), format(tot, ",")))
    if not te_features:
        sys.stderr.write(
            "te_matrix: WARNING - no TE rows found. TEcount emits TE features as\n"
            "  'name:family:class'; none were present, which usually means the --TE GTF\n"
            "  was empty or its chromosome names do not match the alignment.\n")
    return len(gene_features), len(te_features)


def _selftest():
    import tempfile
    d = tempfile.mkdtemp()
    a = os.path.join(d, "a.cntTable")
    b = os.path.join(d, "b.cntTable")
    open(a, "w").write("gene/TE\tA\nENSMUSG001\t100\nENSMUSG002\t50\n"
                       "L1Md_T:L1:LINE\t400\nIAPEz-int:ERVK:LTR\t200\n")
    # sample B is missing ENSMUSG002 entirely, which must become 0 rather than dropping
    # the feature or shifting the columns.
    open(b, "w").write("gene/TE\tB\nENSMUSG001\t80\n"
                       "L1Md_T:L1:LINE\t50\nIAPEz-int:ERVK:LTR\t25\n")

    out = [os.path.join(d, x) for x in
           ("combined.tsv", "te.tsv", "gene.tsv", "summary.tsv")]
    n_gene, n_te = run([a, b], ["A", "B"], *out)

    combined = open(out[0]).read()
    te = open(out[1]).read()
    gene = open(out[2]).read()
    summary = open(out[3]).read()

    checks = [
        ("gene and TE rows separated", n_gene == 2 and n_te == 2),
        ("TE matrix holds only TE rows", "ENSMUSG001" not in te and "L1Md_T" in te),
        ("gene matrix holds only genes", "L1Md_T" not in gene and "ENSMUSG001" in gene),
        ("missing feature filled with 0", "ENSMUSG002\t50\t0" in combined),
        ("class rollup present", "class\tLINE" in summary and "class\tLTR" in summary),
        ("family rollup present", "family\tL1" in summary and "family\tERVK" in summary),
        # A 600/750 TE share in sample A: 400+200 TE against 150 gene counts.
        ("TE share computed", "0.80000" in summary),
    ]

    # Negative control: a gene ID containing one colon must not be read as a TE.
    checks.append(("single-colon ID is not a TE", not is_te("weird:gene")))
    checks.append(("three-part ID is a TE", is_te("L1Md_T:L1:LINE")))

    ok = True
    for label, passed in checks:
        print("  %s  %s" % ("PASS" if passed else "FAIL", label))
        ok = ok and passed
    print("\n%s" % ("SELFTEST PASSED" if ok else "SELFTEST FAILED"))
    return 0 if ok else 1


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        sys.exit(_selftest())
    try:
        snakemake  # noqa: F821
    except NameError:
        print(__doc__)
        sys.exit("run via snakemake, or with --selftest")
    run(list(snakemake.input.counts),            # noqa: F821
        list(snakemake.params.samples),          # noqa: F821
        snakemake.output.combined,               # noqa: F821
        snakemake.output.te_only,                # noqa: F821
        snakemake.output.gene_only,              # noqa: F821
        snakemake.output.summary)                # noqa: F821
