#!/usr/bin/env python3
"""Split single-end FASTQ into simulated paired-end by halving each read and
reverse-complementing the second half. Sequence-lossless aside from a small
configurable gap and optional per-read jitter; no synthetic bases or simulated
errors are introduced."""

import argparse
import gzip
import random
import sys
from pathlib import Path

# IUPAC complement table including ambiguity codes. Lowercase preserved so
# soft-masked bases stay soft-masked after revcomp.
COMPLEMENT = str.maketrans("ACGTNacgtnRYSWKMBDHVryswkmbdhv",
                           "TGCANtgcanYRSWMKVHDByrswmkvhdb")


def revcomp(seq: str) -> str:
    # translate() complements; [::-1] reverses. Combined = reverse complement.
    return seq.translate(COMPLEMENT)[::-1]


def smart_open(path: str, mode: str):
    # Transparent .gz handling plus '-' for stdin/stdout so the script can be
    # piped through pigz or similar by the caller if higher throughput is needed.
    if path == "-":
        return sys.stdin if "r" in mode else sys.stdout
    if path.endswith(".gz"):
        return gzip.open(path, mode, compresslevel=4)
    return open(path, mode)


def rename_for_pair(header: str, mate: int) -> str:
    """Rewrite a FASTQ header so R1 and R2 records are recognised as a pair.

    Handles three common header conventions:
      1. Illumina with comment: `@<name> 1:N:0:INDEX`  -> flip the read-number
         field in the comment.
      2. Old-style suffix:      `@<name>/1`            -> swap /1 <-> /2.
      3. Bare header:           `@<name>`              -> append /1 or /2.
    """
    header = header.rstrip("\n")
    if not header.startswith("@"):
        raise ValueError(f"Malformed FASTQ header: {header!r}")
    if " " in header:
        # Illumina-style: split off the colon-delimited comment after the name.
        base, rest = header.split(" ", 1)
        parts = rest.split(":")
        # If the first comment field is already '1' or '2', overwrite it.
        # Otherwise prepend the mate number so aligners can still pair on name.
        if parts and parts[0] in ("1", "2"):
            parts[0] = str(mate)
            return f"{base} {':'.join(parts)}"
        return f"{base} {mate}:{rest}"
    if header.endswith(("/1", "/2")):
        return header[:-2] + f"/{mate}"
    return f"{header}/{mate}"


def split_record(seq: str, qual: str, gap_target: int, jitter_target: int,
                 min_mate_len: int, rng: random.Random):
    """Split one SE record into R1 and reverse-complemented R2.

    Adaptive: gap and jitter shrink per-read so short reads (e.g. heavily
    quality-trimmed) still produce a valid pair as long as the read can hold
    two mates of `min_mate_len`. Returns None when the read is too short.
    """
    L = len(seq)
    # max_extra = bases available beyond the two minimum-length mates.
    # Anything we want to drop (gap, jitter) has to fit in here.
    max_extra = L - 2 * min_mate_len
    if max_extra < 0:
        return None

    # Use the configured gap when the read can afford it; shrink to fit otherwise.
    gap = min(gap_target, max_extra)
    # Jitter fills whatever budget remains after the gap, up to the configured max.
    # Per-read random draw is what makes R2's reference-end position vary across
    # otherwise-identical SE reads, which is what prevents MarkDuplicates from
    # collapsing them into one cluster.
    j_max = min(jitter_target, max_extra - gap)
    j = rng.randint(0, j_max) if j_max > 0 else 0
    drop = gap + j

    # Remaining bases are split evenly between the two mates; odd byte to R1.
    usable = L - drop
    r1_len = (usable + 1) // 2
    r2_len = usable - r1_len

    # R1 is the literal 5' prefix of the SE read - sequence and quality unchanged.
    r1_seq = seq[:r1_len]
    r1_qual = qual[:r1_len]

    # R2 starts after R1 + dropped middle (gap + jitter). It is reverse-complemented
    # so the pair aligns on opposite strands; quality is REVERSED but NOT
    # complemented (quality scores are not strand-specific).
    r2_start = r1_len + drop
    r2_seq_fwd = seq[r2_start:r2_start + r2_len]
    r2_qual_fwd = qual[r2_start:r2_start + r2_len]
    r2_seq = revcomp(r2_seq_fwd)
    r2_qual = r2_qual_fwd[::-1]
    return r1_seq, r1_qual, r2_seq, r2_qual


def main():
    ap = argparse.ArgumentParser(
        description="Lossy-symmetric split of single-end FASTQ into simulated "
                    "paired-end by halving each read and reverse-complementing R2."
    )
    ap.add_argument("--in", dest="infile", default="-",
                    help="Input SE FASTQ (.fastq or .fastq.gz). '-' for stdin.")
    ap.add_argument("--r1", required=True, help="Output R1 FASTQ(.gz)")
    ap.add_argument("--r2", required=True, help="Output R2 FASTQ(.gz)")
    ap.add_argument("--gap-size", type=int, default=10,
                    help="Target gap between R1 and R2. Auto-shrinks per read "
                         "down to 0 if the read is too short to hold both "
                         "mates plus the full gap.")
    ap.add_argument("--jitter-max", type=int, default=0,
                    help="Target max random extra bp dropped at R2 start per "
                         "read (drawn uniformly from [0, jitter_max]). "
                         "Auto-shrinks if the read can't fit both mates plus "
                         "gap plus jitter. Breaks MarkDuplicates grouping "
                         "when >0.")
    ap.add_argument("--seed", type=int, default=42,
                    help="RNG seed for jitter (reproducibility).")
    ap.add_argument("--min-mate-len", type=int, default=30,
                    help="Minimum length per mate; shorter records are dropped.")
    ap.add_argument("--stats", default=None,
                    help="Optional path to write summary stats (TSV).")
    args = ap.parse_args()

    if args.jitter_max < 0:
        sys.exit("--jitter-max must be >= 0")
    # Seeded RNG = deterministic output for the same input. Important: re-running
    # the pipeline must produce byte-identical FASTQs so downstream variant calls
    # are reproducible.
    rng = random.Random(args.seed)

    # Running counters for the optional --stats TSV and the stderr summary line.
    n_in = 0
    n_out = 0
    n_dropped = 0
    total_bp_in = 0
    total_bp_out = 0

    with smart_open(args.infile, "rt") as fin, \
            smart_open(args.r1, "wt") as f1, \
            smart_open(args.r2, "wt") as f2:
        # Stream record-by-record; no whole-file buffering, so memory is O(1)
        # regardless of input size.
        while True:
            h = fin.readline()
            if not h:
                break  # EOF
            s = fin.readline().rstrip("\n")
            p = fin.readline()
            q = fin.readline().rstrip("\n")
            # Cheap structural validation - cheap enough to do every record.
            if not q or not p.startswith("+"):
                sys.exit("FASTQ format error near record "
                         f"{n_in + 1}: missing '+' or quality line")
            n_in += 1
            total_bp_in += len(s)
            if len(s) != len(q):
                sys.exit(f"Seq/qual length mismatch at record {n_in}")

            split = split_record(s, q, args.gap_size, args.jitter_max,
                                 args.min_mate_len, rng)
            if split is None:
                # Read too short to produce two mates of min_mate_len.
                n_dropped += 1
                continue
            r1_seq, r1_qual, r2_seq, r2_qual = split
            # Rename per mate so name-based pair detection in aligners works.
            h1 = rename_for_pair(h, 1)
            h2 = rename_for_pair(h, 2)
            f1.write(f"{h1}\n{r1_seq}\n+\n{r1_qual}\n")
            f2.write(f"{h2}\n{r2_seq}\n+\n{r2_qual}\n")
            n_out += 1
            total_bp_out += len(r1_seq) + len(r2_seq)

    # Stats TSV is consumed by MultiQC custom-content or by humans; one row per run.
    if args.stats:
        Path(args.stats).write_text(
            "reads_in\treads_out\treads_dropped\tbp_in\tbp_out\tbp_retained_pct\n"
            f"{n_in}\t{n_out}\t{n_dropped}\t{total_bp_in}\t{total_bp_out}\t"
            f"{(100 * total_bp_out / total_bp_in) if total_bp_in else 0:.2f}\n"
        )

    # Compact one-liner on stderr for Nextflow's .command.err log.
    sys.stderr.write(
        f"[split_se_to_pe] in={n_in} out_pairs={n_out} dropped={n_dropped} "
        f"bp_in={total_bp_in} bp_out={total_bp_out}\n"
    )


if __name__ == "__main__":
    main()
