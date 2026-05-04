#!/usr/bin/env python3
import argparse
import gzip
import random
import sys
from pathlib import Path

COMPLEMENT = str.maketrans("ACGTNacgtnRYSWKMBDHVryswkmbdhv",
                           "TGCANtgcanYRSWMKVHDByrswmkvhdb")


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def smart_open(path: str, mode: str):
    if path == "-":
        return sys.stdin if "r" in mode else sys.stdout
    if path.endswith(".gz"):
        return gzip.open(path, mode, compresslevel=4)
    return open(path, mode)


def rename_for_pair(header: str, mate: int) -> str:
    header = header.rstrip("\n")
    if not header.startswith("@"):
        raise ValueError(f"Malformed FASTQ header: {header!r}")
    if " " in header:
        base, rest = header.split(" ", 1)
        parts = rest.split(":")
        if parts and parts[0] in ("1", "2"):
            parts[0] = str(mate)
            return f"{base} {':'.join(parts)}"
        return f"{base} {mate}:{rest}"
    if header.endswith(("/1", "/2")):
        return header[:-2] + f"/{mate}"
    return f"{header}/{mate}"


def split_record(seq: str, qual: str, gap_target: int, jitter_target: int,
                 min_mate_len: int, rng: random.Random):
    L = len(seq)
    max_extra = L - 2 * min_mate_len
    if max_extra < 0:
        return None
    gap = min(gap_target, max_extra)
    j_max = min(jitter_target, max_extra - gap)
    j = rng.randint(0, j_max) if j_max > 0 else 0
    drop = gap + j
    usable = L - drop
    r1_len = (usable + 1) // 2
    r2_len = usable - r1_len
    r1_seq = seq[:r1_len]
    r1_qual = qual[:r1_len]
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
    rng = random.Random(args.seed)

    n_in = 0
    n_out = 0
    n_dropped = 0
    total_bp_in = 0
    total_bp_out = 0

    with smart_open(args.infile, "rt") as fin, \
            smart_open(args.r1, "wt") as f1, \
            smart_open(args.r2, "wt") as f2:
        while True:
            h = fin.readline()
            if not h:
                break
            s = fin.readline().rstrip("\n")
            p = fin.readline()
            q = fin.readline().rstrip("\n")
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
                n_dropped += 1
                continue
            r1_seq, r1_qual, r2_seq, r2_qual = split
            h1 = rename_for_pair(h, 1)
            h2 = rename_for_pair(h, 2)
            f1.write(f"{h1}\n{r1_seq}\n+\n{r1_qual}\n")
            f2.write(f"{h2}\n{r2_seq}\n+\n{r2_qual}\n")
            n_out += 1
            total_bp_out += len(r1_seq) + len(r2_seq)

    if args.stats:
        Path(args.stats).write_text(
            "reads_in\treads_out\treads_dropped\tbp_in\tbp_out\tbp_retained_pct\n"
            f"{n_in}\t{n_out}\t{n_dropped}\t{total_bp_in}\t{total_bp_out}\t"
            f"{(100 * total_bp_out / total_bp_in) if total_bp_in else 0:.2f}\n"
        )

    sys.stderr.write(
        f"[split_se_to_pe] in={n_in} out_pairs={n_out} dropped={n_dropped} "
        f"bp_in={total_bp_in} bp_out={total_bp_out}\n"
    )


if __name__ == "__main__":
    main()
