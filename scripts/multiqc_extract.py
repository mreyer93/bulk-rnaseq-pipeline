"""Recover per-sample metrics from a MultiQC HTML report.

Sequencing cores routinely hand over MultiQC HTML and nothing else. The numbers are in
there, but they are LZ-string-compressed inside a script tag, so they cannot be grepped
and are not obviously reusable. This turns a report back into tab-separated tables you can
actually analyse: one file per plot, samples as rows, series as columns.

Works with MultiQC 1.x reports that embed `mqc_compressed_plotdata` (1.9 through at least
1.14). Older reports that store plain JSON in `mqc_plots` are handled too.

Run standalone:
    python scripts/multiqc_extract.py report.html --outdir metrics/
    python scripts/multiqc_extract.py report.html --list
    python scripts/multiqc_extract.py --selftest

Needs `lzstring` (pip install lzstring) only for the compressed flavour.
"""

# NB: no `from __future__ import annotations` - see tx2gene.py for why.

import argparse
import json
import os
import re
import sys

_COMPRESSED = re.compile(r'id="mqc_compressed_plotdata"[^>]*>(.*?)</script>', re.S)
_PLAIN = re.compile(r'mqc_plots\s*=\s*(\{.*?\});\s*\n', re.S)


def _decompress(payload):
    try:
        import lzstring
    except ImportError:
        raise SystemExit(
            "This report stores its data LZ-string-compressed; install the decoder with:\n"
            "    pip install lzstring")
    return lzstring.LZString().decompressFromBase64(payload.strip())


def load_plots(html_path):
    """Return the MultiQC plot-data dict, whichever way the report stored it."""
    with open(html_path, encoding="utf-8", errors="replace") as fh:
        html = fh.read()

    m = _COMPRESSED.search(html)
    if m:
        text = _decompress(m.group(1))
        if not text:
            raise SystemExit("found a compressed payload but it did not decode")
        return json.loads(text)

    m = _PLAIN.search(html)
    if m:
        try:
            return json.loads(m.group(1))
        except ValueError:
            pass
    raise SystemExit("no MultiQC plot data found in {}".format(html_path))


def bar_table(plot):
    """Flatten a MultiQC bar_graph into (sample_names, {series: {sample: value}})."""
    samples = plot.get("samples") or []
    names = samples[0] if samples and isinstance(samples[0], list) else samples
    datasets = plot.get("datasets") or []
    if not datasets:
        return names, {}
    first = datasets[0]
    series = {}
    if isinstance(first, list):
        for s in first:
            if isinstance(s, dict) and "name" in s:
                series[s["name"]] = dict(zip(names, s.get("data", [])))
    return names, series


def write_tsv(names, series, out_path, add_total=True):
    cols = list(series)
    with open(out_path, "w") as out:
        header = ["sample"] + cols + (["total"] if add_total else [])
        out.write("\t".join(header) + "\n")
        for s in names:
            vals = [series[c].get(s, "") for c in cols]
            row = [s] + [("" if v == "" else repr(v) if isinstance(v, float) else str(v))
                         for v in vals]
            if add_total:
                nums = [v for v in vals if isinstance(v, (int, float))]
                row.append(str(sum(nums)) if nums else "")
            out.write("\t".join(row) + "\n")
    return len(names), len(cols)


def extract(html_path, outdir, quiet=False):
    plots = load_plots(html_path)
    if outdir and not os.path.isdir(outdir):
        os.makedirs(outdir)
    written = []
    for key, plot in plots.items():
        if plot.get("plot_type") != "bar_graph":
            continue
        names, series = bar_table(plot)
        if not names or not series:
            continue
        out_path = os.path.join(outdir or ".", "{}.tsv".format(key))
        n, c = write_tsv(names, series, out_path)
        written.append((out_path, n, c))
        if not quiet:
            sys.stderr.write("wrote {}  ({} samples x {} series)\n".format(out_path, n, c))
    if not written and not quiet:
        sys.stderr.write("no bar-graph plots found; available plots: {}\n"
                         .format(", ".join(plots)))
    return written


def _selftest():
    """Round-trip a synthetic MultiQC payload through both storage flavours."""
    import io
    plot = {
        "plot_type": "bar_graph",
        "samples": [["sampleA", "sampleB"]],
        "datasets": [[{"name": "Assigned", "data": [80, 60]},
                      {"name": "Unassigned", "data": [20, 40]}]],
    }
    names, series = bar_table(plot)
    checks = [
        ("sample names recovered", names == ["sampleA", "sampleB"]),
        ("series recovered", set(series) == {"Assigned", "Unassigned"}),
        ("values aligned to samples", series["Assigned"]["sampleB"] == 60),
    ]

    buf = io.StringIO()
    cols = list(series)
    buf.write("\t".join(["sample"] + cols + ["total"]) + "\n")
    for s in names:
        vals = [series[c][s] for c in cols]
        buf.write("\t".join([s] + [str(v) for v in vals] + [str(sum(vals))]) + "\n")
    text = buf.getvalue()
    checks.append(("totals computed", "\t100\n" in text))

    # A plain-JSON report must parse through the same path.
    html = 'x<script>mqc_plots = {"p": %s};\n</script>' % json.dumps(plot)
    m = _PLAIN.search(html)
    checks.append(("plain mqc_plots flavour parses",
                   bool(m) and json.loads(m.group(1))["p"]["plot_type"] == "bar_graph"))

    # Negative control: a report with neither flavour must fail loudly, not silently.
    import tempfile
    with tempfile.NamedTemporaryFile("w", suffix=".html", delete=False) as fh:
        fh.write("<html>nothing here</html>")
        empty = fh.name
    try:
        load_plots(empty)
        checks.append(("missing data raises", False))
    except SystemExit:
        checks.append(("missing data raises", True))
    finally:
        os.unlink(empty)

    ok = True
    for label, passed in checks:
        print("  {}  {}".format("PASS" if passed else "FAIL", label))
        ok = ok and passed
    print("\n{}".format("SELFTEST PASSED" if ok else "SELFTEST FAILED"))
    return 0 if ok else 1


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("report", nargs="?", help="MultiQC HTML report")
    p.add_argument("--outdir", default=".", help="directory for the TSV files")
    p.add_argument("--list", action="store_true", help="list plots and exit")
    p.add_argument("--selftest", action="store_true", help="run built-in checks and exit")
    args = p.parse_args(argv)

    if args.selftest:
        return _selftest()
    if not args.report:
        p.error("give a MultiQC HTML report (or --selftest)")
    if args.list:
        for k, v in load_plots(args.report).items():
            n = len((v.get("samples") or [[]])[0]) if v.get("samples") else 0
            print("{}\t{}\t{} samples".format(k, v.get("plot_type"), n))
        return 0
    return 0 if extract(args.report, args.outdir) else 1


if __name__ == "__main__":
    sys.exit(main())
