#!/usr/bin/env python3
import argparse
import csv
import sys
from pathlib import Path

REQUIRED = ["sample", "fastq"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input")
    ap.add_argument("output")
    args = ap.parse_args()

    with open(args.input) as fh:
        reader = csv.DictReader(fh)
        missing = [c for c in REQUIRED if c not in (reader.fieldnames or [])]
        if missing:
            sys.exit(f"Samplesheet missing columns: {missing}")
        rows = list(reader)

    if not rows:
        sys.exit("Samplesheet is empty")

    seen = set()
    for i, row in enumerate(rows, start=2):
        sid = row["sample"].strip()
        fq = row["fastq"].strip()
        if not sid or not fq:
            sys.exit(f"Row {i}: empty sample or fastq")
        if sid in seen:
            sys.exit(f"Row {i}: duplicate sample ID {sid!r}")
        seen.add(sid)
        if not (fq.endswith(".fastq.gz") or fq.endswith(".fq.gz")
                or fq.endswith(".fastq") or fq.endswith(".fq")):
            sys.exit(f"Row {i}: {fq!r} does not look like FASTQ")

    Path(args.output).write_text(
        "sample,fastq\n" +
        "\n".join(f"{r['sample'].strip()},{r['fastq'].strip()}" for r in rows) +
        "\n"
    )


if __name__ == "__main__":
    main()
