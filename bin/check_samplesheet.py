#!/usr/bin/env python3
"""Validate and normalise the input samplesheet.

Accepts either a `reads` column (preferred) or a legacy `fastq` column, and
emits a normalised two-column CSV (`sample,reads`) for the pipeline to consume.
Read files may be FASTQ, BAM, or CRAM.
"""

import argparse
import csv
import sys
from pathlib import Path

# Accepted input types. FASTQ is used as-is; BAM/CRAM are unpacked to FASTQ by
# the CRAM_TO_FASTQ process before splitting.
FASTQ_EXTS = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
ALIGNED_EXTS = (".bam", ".cram")
READ_EXTS = FASTQ_EXTS + ALIGNED_EXTS


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input")
    ap.add_argument("output")
    args = ap.parse_args()

    with open(args.input) as fh:
        reader = csv.DictReader(fh)
        cols = reader.fieldnames or []
        if "sample" not in cols:
            sys.exit("Samplesheet missing required column: 'sample'")
        # 'reads' is the current column name; 'fastq' is kept as an alias so
        # samplesheets written before BAM/CRAM support keep working.
        read_col = "reads" if "reads" in cols else ("fastq" if "fastq" in cols else None)
        if read_col is None:
            sys.exit("Samplesheet missing required column: 'reads' (or legacy 'fastq')")
        rows = list(reader)

    if not rows:
        sys.exit("Samplesheet is empty")

    seen = set()
    out = []
    for i, row in enumerate(rows, start=2):
        sid = (row["sample"] or "").strip()
        path = (row[read_col] or "").strip()
        if not sid or not path:
            sys.exit(f"Row {i}: empty sample or {read_col}")
        if " " in sid:
            sys.exit(f"Row {i}: sample name {sid!r} cannot contain spaces")
        if sid in seen:
            sys.exit(f"Row {i}: duplicate sample ID {sid!r}")
        seen.add(sid)
        if not path.lower().endswith(READ_EXTS):
            sys.exit(
                f"Row {i}: {path!r} is not a recognised read file. "
                f"Expected one of: {', '.join(READ_EXTS)}"
            )
        out.append((sid, path))

    Path(args.output).write_text(
        "sample,reads\n" + "\n".join(f"{s},{p}" for s, p in out) + "\n"
    )


if __name__ == "__main__":
    main()
