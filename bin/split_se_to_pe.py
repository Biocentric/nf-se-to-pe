#!/usr/bin/env python3
"""Split single-end FASTQ into simulated paired-end by halving each read and
reverse-complementing the second half. Sequence-lossless by default; no
synthetic bases or simulated errors are introduced.

Properties of the output, verified by simulation against ground truth:

* TLEN (fragment length) of each output pair equals the ORIGINAL read length
  exactly, because R1 begins at the read's first base and R2 ends at its last.
  The fragment-length distribution therefore reproduces the input read-length
  distribution. It does NOT depend on --gap-size or --jitter-max.
* Duplicate identity is (read start, read end), which is likewise independent
  of --gap-size and --jitter-max. Marking duplicates on the converted pairs
  recovers true PCR duplicates exactly, and is strictly better than single-end
  duplicate marking, which cannot separate reads sharing a start but differing
  in length.

Consequently --gap-size and --jitter-max discard bases without changing any
downstream signal, and both default to 0. They are retained only so existing
command lines keep working.
"""

import argparse
import gzip
import random
import sys
from pathlib import Path

# IUPAC complement table including ambiguity codes. Lowercase preserved so
# soft-masked bases stay soft-masked after revcomp.
COMPLEMENT = str.maketrans("ACGTNacgtnRYSWKMBDHVryswkmbdhv",
                           "TGCANtgcanYRSWMKVHDByrswmkvhdb")

# A split falling strictly inside a microsatellite truncates the repeat tract,
# corrupting repeat-length evidence used by indel callers and by microsatellite
# jitter models (e.g. REDUX feeding SAGE). These thresholds define what counts
# as a tract worth protecting.
MS_MIN_HOMOPOLYMER = 5
MS_MIN_UNITS = 3
MS_MAX_UNIT = 4


def revcomp(seq: str) -> str:
    # translate() complements; [::-1] reverses. Combined = reverse complement.
    return seq.translate(COMPLEMENT)[::-1]


def smart_open(path: str, mode: str):
    # Transparent .gz handling plus '-' for stdin/stdout so the caller can pipe
    # through pigz or similar if higher throughput is needed.
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


def _cuts_homopolymer(seq: str, k: int) -> bool:
    """Does the boundary before index k sit strictly inside a homopolymer run?"""
    n = len(seq)
    if k <= 0 or k >= n or seq[k - 1] != seq[k]:
        return False
    c = seq[k]
    i = k - 1
    while i > 0 and seq[i - 1] == c:
        i -= 1
    j = k
    while j + 1 < n and seq[j + 1] == c:
        j += 1
    return (j - i + 1) >= MS_MIN_HOMOPOLYMER


def _cuts_tandem(seq: str, k: int) -> bool:
    """Does the boundary before index k sit strictly inside a 2-4 bp tandem repeat?

    The candidate unit is anchored at seq[k-u:k]; because a tandem repeat is
    periodic, anchoring at the boundary finds the tract whatever its phase.
    Work is bounded by the tract length, so this stays cheap per read.
    """
    n = len(seq)
    for u in range(2, MS_MAX_UNIT + 1):
        start = k - u
        if start < 0 or k + u > n:
            continue
        unit = seq[start:k]
        i = start
        while i - u >= 0 and seq[i - u:i] == unit:
            i -= u
        j = k
        while j + u <= n and seq[j:j + u] == unit:
            j += u
        # j > k means a further unit lies to the right of the boundary, i.e. the
        # boundary is interior to the tract rather than sitting at its edge.
        if j > k and (j - i) // u >= MS_MIN_UNITS:
            return True
    return False


def cuts_microsatellite(seq: str, k: int) -> bool:
    return _cuts_homopolymer(seq, k) or _cuts_tandem(seq, k)


def choose_split(seq: str, nominal: int, lo: int, hi: int, max_shift: int):
    """Pick a split point, nudging off a microsatellite where cheaply possible.

    Returns (split_point, was_shifted, still_cuts). Shifting the split changes
    only the relative mate lengths; TLEN and duplicate identity are both
    determined by the read's start and end, so the shift costs nothing.
    """
    if max_shift <= 0 or not cuts_microsatellite(seq, nominal):
        return nominal, False, False
    for d in range(1, max_shift + 1):
        for cand in (nominal - d, nominal + d):
            if lo <= cand <= hi and not cuts_microsatellite(seq, cand):
                return cand, True, False
    return nominal, False, True


def split_record(seq: str, qual: str, gap_target: int, jitter_target: int,
                 min_mate_len: int, rng: random.Random, ms_shift: int):
    """Split one SE record into R1 and reverse-complemented R2.

    Adaptive: gap and jitter shrink per-read so short reads (e.g. heavily
    quality-trimmed) still produce a valid pair as long as the read can hold
    two mates of `min_mate_len`. Returns None when the read is too short.
    """
    L = len(seq)
    # max_extra = bases available beyond the two minimum-length mates. Anything
    # we drop (gap, jitter) has to fit in here.
    max_extra = L - 2 * min_mate_len
    if max_extra < 0:
        return None

    # Both default to 0; see the module docstring for why they are no-ops.
    gap = min(gap_target, max_extra)
    j_max = min(jitter_target, max_extra - gap)
    j = rng.randint(0, j_max) if j_max > 0 else 0
    drop = gap + j

    usable = L - drop
    nominal = (usable + 1) // 2
    r1_len, shifted, still_cuts = choose_split(
        seq, nominal, min_mate_len, L - drop - min_mate_len, ms_shift)
    r2_len = usable - r1_len

    # R1 is the literal 5' prefix of the SE read - sequence and quality unchanged.
    r1_seq = seq[:r1_len]
    r1_qual = qual[:r1_len]

    # R2 starts after R1 plus any dropped middle. It is reverse-complemented so
    # the pair aligns on opposite strands; quality is REVERSED but NOT
    # complemented (quality scores are not strand-specific).
    r2_start = r1_len + drop
    r2_seq = revcomp(seq[r2_start:r2_start + r2_len])
    r2_qual = qual[r2_start:r2_start + r2_len][::-1]
    return r1_seq, r1_qual, r2_seq, r2_qual, shifted, still_cuts


def main():
    ap = argparse.ArgumentParser(
        description="Split single-end FASTQ into simulated paired-end by "
                    "halving each read and reverse-complementing R2."
    )
    ap.add_argument("--in", dest="infile", default="-",
                    help="Input SE FASTQ (.fastq or .fastq.gz). '-' for stdin.")
    ap.add_argument("--r1", required=True, help="Output R1 FASTQ(.gz)")
    ap.add_argument("--r2", required=True, help="Output R2 FASTQ(.gz)")
    ap.add_argument("--gap-size", type=int, default=0,
                    help="Bases dropped between the mates. Default 0: the gap "
                         "changes neither fragment length nor duplicate "
                         "identity, so it only discards data. Retained for "
                         "backwards compatibility.")
    ap.add_argument("--jitter-max", type=int, default=0,
                    help="Extra random bases dropped at the R2 start. Default "
                         "0: R2's 5' reference position is the read end "
                         "regardless of jitter, so this cannot affect duplicate "
                         "grouping. Retained for backwards compatibility.")
    ap.add_argument("--ms-shift", type=int, default=6,
                    help="Shift the split point by up to this many bases to "
                         "avoid cutting through a homopolymer or short tandem "
                         "repeat. 0 disables.")
    ap.add_argument("--seed", type=int, default=42,
                    help="RNG seed for jitter (reproducibility).")
    ap.add_argument("--min-mate-len", type=int, default=30,
                    help="Minimum length per mate; shorter records are dropped.")
    ap.add_argument("--stats", default=None,
                    help="Optional path to write summary stats (TSV).")
    ap.add_argument("--mqc", default=None,
                    help="Optional path to write a MultiQC custom-content table.")
    ap.add_argument("--sample", default="sample",
                    help="Sample name used in the MultiQC table.")
    args = ap.parse_args()

    if args.jitter_max < 0:
        sys.exit("--jitter-max must be >= 0")
    if args.gap_size < 0:
        sys.exit("--gap-size must be >= 0")
    # Seeded RNG = deterministic output for the same input, so re-running the
    # pipeline gives byte-identical FASTQs and reproducible variant calls.
    rng = random.Random(args.seed)

    n_in = n_out = n_dropped = 0
    total_bp_in = total_bp_out = 0
    n_shifted = n_residual = 0

    with smart_open(args.infile, "rt") as fin, \
            smart_open(args.r1, "wt") as f1, \
            smart_open(args.r2, "wt") as f2:
        # Stream record-by-record; memory is O(1) regardless of input size.
        while True:
            h = fin.readline()
            if not h:
                break  # EOF
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
                                 args.min_mate_len, rng, args.ms_shift)
            if split is None:
                # Too short to yield two mates of min_mate_len.
                n_dropped += 1
                continue
            r1_seq, r1_qual, r2_seq, r2_qual, shifted, still_cuts = split
            n_shifted += shifted
            n_residual += still_cuts
            f1.write(f"{rename_for_pair(h, 1)}\n{r1_seq}\n+\n{r1_qual}\n")
            f2.write(f"{rename_for_pair(h, 2)}\n{r2_seq}\n+\n{r2_qual}\n")
            n_out += 1
            total_bp_out += len(r1_seq) + len(r2_seq)

    pct = (100 * total_bp_out / total_bp_in) if total_bp_in else 0.0

    if args.stats:
        Path(args.stats).write_text(
            "reads_in\treads_out\treads_dropped\tbp_in\tbp_out\t"
            "bp_retained_pct\tms_splits_shifted\tms_splits_residual\n"
            f"{n_in}\t{n_out}\t{n_dropped}\t{total_bp_in}\t{total_bp_out}\t"
            f"{pct:.2f}\t{n_shifted}\t{n_residual}\n"
        )

    if args.mqc:
        Path(args.mqc).write_text(
            "# id: 'se_to_pe'\n"
            "# section_name: 'SE to PE conversion'\n"
            "# description: 'Summary of the single-end to paired-end split. "
            "Fragment length of the output pairs equals the original read "
            "length, so the converted library inherits the input read-length "
            "distribution.'\n"
            "# plot_type: 'table'\n"
            "# pconfig:\n"
            "#     id: 'se_to_pe_table'\n"
            "#     namespace: 'se-to-pe'\n"
            "Sample\tReads in\tPairs out\tReads dropped\tBases retained (%)\t"
            "MS splits avoided\tMS splits residual\n"
            f"{args.sample}\t{n_in}\t{n_out}\t{n_dropped}\t{pct:.2f}\t"
            f"{n_shifted}\t{n_residual}\n"
        )

    # Compact one-liner on stderr for Nextflow's .command.err log.
    sys.stderr.write(
        f"[split_se_to_pe] in={n_in} out_pairs={n_out} dropped={n_dropped} "
        f"bp_in={total_bp_in} bp_out={total_bp_out} ({pct:.2f}%) "
        f"ms_shifted={n_shifted} ms_residual={n_residual}\n"
    )


if __name__ == "__main__":
    main()
