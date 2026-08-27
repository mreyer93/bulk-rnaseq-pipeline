"""Sample sheet parsing, normalisation and validation for the bulk RNA-seq pipeline.

RNA-seq sample sheets arrive in wildly inconsistent shapes: CSV or TSV, columns called
`sample` / `Sample_ID` / `sampleName`, reads called `fastq_1` / `R1` / `read1` / `forward`,
single- or paired-end, sometimes several FASTQ pairs per sample (one per lane). This module
accepts all of those, normalises them to one internal schema, and fails with a specific,
actionable message when it can't.

Usable two ways:
  * imported by the workflow (`load_samplesheet`)
  * run directly as a pre-flight check:  python scripts/samplesheet.py samples.csv
"""

from __future__ import annotations

import csv
import os
import re
import sys
from collections import defaultdict, OrderedDict

# Column aliases, lowercased and stripped of non-alphanumerics before matching, so
# "Sample ID", "sample_id" and "SampleID" all collapse to the same key.
_ALIASES = {
    "sample":    ["sample", "sampleid", "samplename", "name", "id", "sampleaccession"],
    "fastq_1":   ["fastq1", "fastqr1", "r1", "read1", "reads1", "forward", "fq1",
                  "file1", "fastq", "reads", "filename1"],
    "fastq_2":   ["fastq2", "fastqr2", "r2", "read2", "reads2", "reverse", "fq2",
                  "file2", "filename2"],
    "condition": ["condition", "group", "treatment", "class", "phenotype", "status"],
    "batch":     ["batch", "run", "lane", "flowcell", "plate"],
    "strandedness": ["strandedness", "strand", "library_strandedness", "librarytype"],
}

_VALID_STRANDEDNESS = {"auto", "unstranded", "forward", "reverse"}

# A sample name that survives into filenames, R factor levels and plot labels
_SAFE_SAMPLE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


class SampleSheetError(ValueError):
    """Raised with a message intended to be read directly by the user."""


def _norm_key(k: str) -> str:
    return re.sub(r"[^a-z0-9]", "", (k or "").strip().lower())


def _sniff_delimiter(path: str) -> str:
    with open(path, newline="", encoding="utf-8-sig") as fh:
        head = fh.readline()
    if not head:
        raise SampleSheetError(f"{path} is empty.")
    # Prefer whichever of tab/comma/semicolon appears most in the header line
    counts = {d: head.count(d) for d in ("\t", ",", ";")}
    delim = max(counts, key=counts.get)
    if counts[delim] == 0:
        raise SampleSheetError(
            f"Could not find a tab, comma or semicolon separator in the header of {path}.\n"
            f"Header was: {head.strip()!r}"
        )
    return delim


def _map_columns(fieldnames):
    """Map the sheet's actual column names onto our internal schema."""
    mapping = {}
    unmapped = []
    for raw in fieldnames:
        key = _norm_key(raw)
        for canonical, aliases in _ALIASES.items():
            if key == _norm_key(canonical) or key in aliases:
                # first column to claim a canonical name wins; later ones stay extra
                if canonical not in mapping.values():
                    mapping[raw] = canonical
                    break
        else:
            unmapped.append(raw)
    return mapping, unmapped


def load_samplesheet(path: str, check_files: bool = True, base_dir: str | None = None):
    """Parse and validate a sample sheet.

    Returns (samples, extra_columns) where `samples` is an OrderedDict keyed by sample
    name; each value is a dict with keys: name, units (list of {fastq_1, fastq_2}),
    paired (bool), strandedness, and any extra metadata columns (condition, batch, ...).
    """
    if not os.path.exists(path):
        raise SampleSheetError(f"Sample sheet not found: {path}")

    delim = _sniff_delimiter(path)
    with open(path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        if not reader.fieldnames:
            raise SampleSheetError(f"{path} has no header row.")
        mapping, unmapped = _map_columns(reader.fieldnames)
        canonical = set(mapping.values())

        if "sample" not in canonical:
            raise SampleSheetError(
                f"No sample-name column found in {path}.\n"
                f"Columns present: {reader.fieldnames}\n"
                f"Name one of them 'sample' (or sample_id / sample_name / name / id)."
            )
        if "fastq_1" not in canonical:
            raise SampleSheetError(
                f"No FASTQ column found in {path}.\n"
                f"Columns present: {reader.fieldnames}\n"
                f"Name one of them 'fastq_1' (or R1 / read1 / forward / fastq)."
            )

        rows = list(reader)

    if not rows:
        raise SampleSheetError(f"{path} has a header but no data rows.")

    base_dir = base_dir or os.path.dirname(os.path.abspath(path))

    def resolve(p):
        p = (p or "").strip()
        if not p:
            return ""
        return p if os.path.isabs(p) else os.path.normpath(os.path.join(base_dir, p))

    samples = OrderedDict()
    errors = []
    # `extra` = user metadata columns we pass through to the design matrix
    extra_cols = [c for c in mapping.values()
                  if c not in ("sample", "fastq_1", "fastq_2", "strandedness")]
    extra_cols += [c for c in unmapped]

    for i, row in enumerate(rows, start=2):  # start=2 -> header is line 1
        rec = {}
        for raw, canon in mapping.items():
            rec[canon] = (row.get(raw) or "").strip()
        for raw in unmapped:
            rec[raw] = (row.get(raw) or "").strip()

        name = rec.get("sample", "")
        if not name:
            errors.append(f"line {i}: empty sample name")
            continue
        if not _SAFE_SAMPLE.match(name):
            errors.append(
                f"line {i}: sample name {name!r} contains characters that break downstream "
                f"filenames/labels; use letters, digits, dot, dash, underscore and start "
                f"with an alphanumeric"
            )
            continue

        f1 = resolve(rec.get("fastq_1", ""))
        f2 = resolve(rec.get("fastq_2", ""))
        if not f1:
            errors.append(f"line {i} ({name}): no FASTQ file given")
            continue
        if check_files:
            for f in filter(None, (f1, f2)):
                if not os.path.exists(f):
                    errors.append(f"line {i} ({name}): FASTQ not found: {f}")

        strand = (rec.get("strandedness") or "auto").strip().lower() or "auto"
        if strand not in _VALID_STRANDEDNESS:
            errors.append(
                f"line {i} ({name}): strandedness {strand!r} is not one of "
                f"{sorted(_VALID_STRANDEDNESS)}"
            )

        if name not in samples:
            samples[name] = {
                "name": name,
                "units": [],
                "paired": bool(f2),
                "strandedness": strand,
            }
            for c in extra_cols:
                samples[name][c] = rec.get(c, "")
        else:
            # repeated sample name = extra sequencing unit (lane/run) to be merged
            if bool(f2) != samples[name]["paired"]:
                errors.append(
                    f"line {i} ({name}): mixes single-end and paired-end rows for the same "
                    f"sample; split them into separate samples or supply both mates"
                )
            for c in extra_cols:
                prev, cur = samples[name].get(c, ""), rec.get(c, "")
                if prev != cur and cur:
                    errors.append(
                        f"line {i} ({name}): conflicting {c!r} value ({prev!r} vs {cur!r}) "
                        f"across rows for the same sample"
                    )

        samples[name]["units"].append({"fastq_1": f1, "fastq_2": f2})

    if errors:
        raise SampleSheetError(
            "Sample sheet problems in {}:\n  - {}".format(path, "\n  - ".join(errors))
        )

    # Duplicate FASTQ paths across different samples is almost always a copy-paste bug
    seen = defaultdict(list)
    for s in samples.values():
        for u in s["units"]:
            for f in filter(None, (u["fastq_1"], u["fastq_2"])):
                seen[f].append(s["name"])
    dupes = {f: sorted(set(n)) for f, n in seen.items() if len(set(n)) > 1}
    if dupes:
        msg = "\n  - ".join(f"{f} used by samples {n}" for f, n in dupes.items())
        raise SampleSheetError(f"The same FASTQ is assigned to multiple samples:\n  - {msg}")

    return samples, extra_cols


def summarise(samples, extra_cols) -> str:
    n_units = sum(len(s["units"]) for s in samples.values())
    pe = sum(1 for s in samples.values() if s["paired"])
    lines = [
        f"{len(samples)} samples, {n_units} FASTQ unit(s)",
        f"  paired-end: {pe}, single-end: {len(samples) - pe}",
    ]
    multi = [s["name"] for s in samples.values() if len(s["units"]) > 1]
    if multi:
        lines.append(f"  merged across multiple units: {', '.join(multi)}")
    if extra_cols:
        lines.append(f"  metadata columns: {', '.join(extra_cols)}")
        for c in extra_cols:
            vals = sorted({s.get(c, "") for s in samples.values() if s.get(c, "")})
            if 0 < len(vals) <= 12:
                lines.append(f"    {c}: {', '.join(vals)}")
    return "\n".join(lines)


def validate_design(samples, extra_cols, design: str, contrast=None):
    """Check a DESeq2 design formula against the sample sheet before the run starts,
    so design mistakes surface in seconds instead of after quantification."""
    problems = []
    terms = [t.strip() for t in re.split(r"[+~*:]", design or "") if t.strip()]
    for t in terms:
        if t not in extra_cols:
            problems.append(
                f"design term {t!r} is not a column in the sample sheet "
                f"(available: {', '.join(extra_cols) or 'none'})"
            )
            continue
        vals = [s.get(t, "") for s in samples.values()]
        if any(v == "" for v in vals):
            problems.append(f"design term {t!r} has blank values for some samples")
        levels = sorted(set(v for v in vals if v))
        if len(levels) < 2:
            problems.append(
                f"design term {t!r} has only one level ({levels}); "
                f"DESeq2 needs at least two groups to compare"
            )
        else:
            counts = {lv: sum(1 for v in vals if v == lv) for lv in levels}
            singles = [lv for lv, n in counts.items() if n < 2]
            if singles:
                problems.append(
                    f"design term {t!r} has group(s) with no replication: "
                    f"{singles} (counts: {counts}). DESeq2 cannot estimate dispersion "
                    f"without replicates"
                )
    if contrast:
        if len(contrast) != 3:
            problems.append(
                f"contrast must be [factor, level_numerator, level_denominator]; got {contrast!r}"
            )
        else:
            factor, num, den = contrast
            if factor not in extra_cols:
                problems.append(f"contrast factor {factor!r} is not a sample sheet column")
            else:
                levels = sorted({s.get(factor, "") for s in samples.values()})
                for lv in (num, den):
                    if lv not in levels:
                        problems.append(
                            f"contrast level {lv!r} not found in column {factor!r} "
                            f"(levels present: {levels})"
                        )
    return problems


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        print("usage: python scripts/samplesheet.py <samplesheet> [design] "
              "[factor,numerator,denominator]")
        return 2
    path = argv[1]
    design = argv[2] if len(argv) > 2 else None
    contrast = argv[3].split(",") if len(argv) > 3 else None
    try:
        samples, extra = load_samplesheet(path, check_files=True)
    except SampleSheetError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    print(summarise(samples, extra))
    if design:
        problems = validate_design(samples, extra, design, contrast)
        if problems:
            print("\nDesign problems:", file=sys.stderr)
            for p in problems:
                print(f"  - {p}", file=sys.stderr)
            return 1
        print(f"\nDesign ~{design} is valid for this sample sheet.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
