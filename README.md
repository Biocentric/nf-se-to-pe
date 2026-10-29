# ![nf-se-to-pe](docs/images/nf-se-to-pe_logo_light.png)

[![GitHub Actions CI Status](https://github.com/Biocentric/nf-se-to-pe/actions/workflows/ci.yml/badge.svg)](https://github.com/Biocentric/nf-se-to-pe/actions/workflows/ci.yml)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed.svg?labelColor=000000)](https://www.docker.com/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049.svg?labelColor=000000)](https://docs.conda.io/en/latest/)

## Introduction

**Biocentric/nf-se-to-pe** converts single-end (SE) sequencing data into
simulated paired-end (PE) FASTQ data, so that downstream tools requiring PE
input (e.g. somatic variant callers that expect `_R1`/`_R2`) can accept it.

Input may be **FASTQ, BAM, or CRAM**. Aligned inputs are unpacked back to
FASTQ first; a CRAM usually decodes **without you supplying the reference
genome it was aligned to** (see [Aligned (BAM/CRAM) input](#aligned-bamcram-input)).

The conversion is **sequence-lossless**: every output base comes from the
original SE read, and no synthetic sequence or error noise is introduced.
Each read is split near its midpoint, and R2 is the reverse-complement of the
second half.

## What changed in 1.2.0

Earlier versions defaulted to `gap_size=10` and `jitter_max=3`, on the
reasoning that the pairs would otherwise have a degenerate fragment-length
distribution and that every pair would be flagged as a PCR duplicate.
**Both of those premises were wrong.** Simulation against ground truth shows:

| Claim (pre-1.2.0) | Measured reality |
| ----------------- | ---------------- |
| "Insert size is a narrow spike at `gap_size`" | **TLEN equals the original read length exactly, in 100% of pairs.** R1 begins at the read's first base and R2 ends at its last, so the fragment-length distribution simply reproduces the input read-length distribution. `gap_size` and `jitter_max` do not affect it at all. |
| "Without jitter, MarkDuplicates flags everything as duplicate" | **Duplicate marking is exact with jitter disabled** — 1895/1895 true duplicates recovered, 0 false positives. Duplicate identity is (read start, read end); R2's 5′ reference position is the read end *regardless* of jitter, so jitter cannot influence grouping. |
| "Halving reads will hurt mappability" | **No measurable penalty.** On a repeat-containing reference, 93.0% of converted mates reached MAPQ ≥ 30 versus 92.5% of the original SE reads — the aligner's mate rescue offsets the shorter length. |

So gap and jitter were discarding data for no benefit. Both now default to
**0**, and the conversion retains 100% of the bases of every read it keeps.
On a real trimmed sample, base retention rose from **86.22% to 94.93%** (the
remaining loss is reads below `2 × min_mate_len`, which cannot be split).

Also new: **`--ms_shift`**, which nudges the split point off microsatellites
(below), and a MultiQC table summarising the conversion.

`gap_size` and `jitter_max` are retained so existing command lines keep
working, but there is no known reason to set either above 0.

## Pipeline summary

1. Samplesheet validation ([`check_samplesheet.py`](bin/check_samplesheet.py))
2. BAM/CRAM → FASTQ unpacking, for aligned inputs only ([`samtools fastq`](https://www.htslib.org/doc/samtools-fasta.html))
3. Input FASTQ QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
4. SE → PE conversion ([`split_se_to_pe.py`](bin/split_se_to_pe.py))
5. Converted FASTQ QC (`FastQC`)
6. Aggregate QC report ([`MultiQC`](https://multiqc.info/))

## Quick start

```bash
nextflow run Biocentric/nf-se-to-pe -profile test,singularity --outdir results_test
```

```bash
nextflow run Biocentric/nf-se-to-pe \
    -profile singularity \
    --input samplesheet.csv \
    --outdir results
```

The defaults are the recommended settings; you should not normally need to
pass any conversion parameters.

## Samplesheet input

FASTQ, BAM, and CRAM may be mixed freely in one samplesheet:

```csv
sample,reads
tumor_01,/path/to/tumor_01.fastq.gz
normal_01,/path/to/normal_01.cram
archive_07,/path/to/archive_07.bam
```

| Column   | Description                                                                                       |
| -------- | ------------------------------------------------------------------------------------------------- |
| `sample` | Custom sample name. Must be unique. Spaces are not allowed.                                       |
| `reads`  | Full path or URL to a single-end read file: `.fastq`/`.fq` (optionally `.gz`), `.bam`, or `.cram`. |

The legacy column name `fastq` is still accepted as an alias for `reads`.

## Aligned (BAM/CRAM) input

BAM and CRAM inputs are converted back to FASTQ with `samtools fastq` before
splitting. Two details are handled for you:

- **Reverse-strand reads are restored to original orientation.** An aligner
  stores reverse-strand reads reverse-complemented with the quality string
  reversed; `samtools fastq` undoes this, so the recovered FASTQ matches the
  pre-alignment reads base-for-base.
- **Secondary (`0x100`) and supplementary (`0x800`) alignments are dropped**
  (`-F 0x900`), since they re-present sequence already carried by the primary
  record.

### CRAM without the reference genome

CRAM is reference-compressed, so decoding normally needs the exact reference
the data was aligned to. You usually do **not** need to supply it. Each `@SQ`
header line carries an `M5` checksum of its reference sequence, and htslib
resolves the reference by that checksum, checking:

1. `--fasta`, if you passed one
2. `REF_CACHE` — a local directory of reference chunks (`--ref_cache`)
3. `REF_PATH` — by default the EBI/ENA reference registry (`--ref_path`)

So a CRAM aligned to any standard genome (GRCh37, GRCh38, hs37d5, …) decodes
over the network with no local FASTA. There is no need to convert
CRAM → BAM → FASTQ yourself first.

- **Compute nodes with internet access**: nothing to do; the default works.
- **Offline / air-gapped clusters**: pass `--fasta`, or pre-populate a cache
  on a connected machine and pass `--ref_cache /shared/hts-ref`.
- **Many samples sharing one reference**: set `--ref_cache` regardless, so
  only the first sample hits the network.
- **Non-standard or in-house reference**: `--fasta` is required.

Paired-end BAM/CRAM input is rejected with an explicit error rather than
being silently mangled.

## Parameters

| Parameter          | Description                                                                                              | Default |
| ------------------ | -------------------------------------------------------------------------------------------------------- | ------- |
| `--input`          | Path to the input samplesheet CSV. **Required.**                                                         | `null`  |
| `--outdir`         | Output directory.                                                                                        | `./results` |
| `--min_mate_len`   | Minimum length per output mate. Reads below `2 × min_mate_len` are dropped.                              | `30`    |
| `--ms_shift`       | Shift the split point by up to this many bases to avoid cutting a homopolymer or short tandem repeat. `0` disables. | `6`     |
| `--gap_size`       | Bases dropped between the mates. **Legacy — has no effect on fragment length or duplicate marking, so it only discards data.** | `0` |
| `--jitter_max`     | Extra random bases dropped at the R2 start. **Legacy — cannot affect duplicate grouping.**               | `0`     |
| `--seed`           | RNG seed for jitter (reproducibility).                                                                   | `42`    |
| `--fasta`          | Reference FASTA for decoding CRAM input. Optional.                                                       | `null`  |
| `--ref_cache`      | Directory caching CRAM reference sequences by MD5.                                                       | `null`  |
| `--ref_path`       | Where htslib looks up CRAM references by MD5.                                                            | EBI/ENA registry |
| `--save_se_fastq`  | Publish the FASTQ unpacked from BAM/CRAM.                                                                | `false` |
| `--singularity_cache_dir` | Override the Singularity image cache directory.                                                   | `~/.singularity-cache/se-to-pe` |

## How the conversion works

For an input read of length `L`:

- Reads shorter than `2 × min_mate_len` are dropped
- `nominal = ⌈L / 2⌉`, then shifted by up to `ms_shift` bases to avoid
  splitting inside a microsatellite
- **R1** = `seq[0 : split]`, original quality
- **R2** = reverse-complement of `seq[split : L]`, quality reversed (not complemented)

so `R1 + revcomp(R2)` reconstructs the original read exactly. Consequences
worth knowing, all verified by simulation:

- **TLEN = L.** The converted library inherits the input read-length
  distribution as its fragment-length distribution.
- **Duplicate identity = (read start, read end).** Marking duplicates on the
  converted pairs is exact, and strictly better than single-end duplicate
  marking, which cannot separate reads that share a start but differ in length.
- **100% properly paired, no singletons.**

Read headers are rewritten to pair cleanly: Illumina-style `1:N:0:…` becomes
`1:N:0:…` / `2:N:0:…`; `/1` suffixes become `/1` / `/2`; bare headers get
`/1` / `/2` appended.

### Microsatellite-aware splitting

A split that lands strictly inside a homopolymer or short tandem repeat
truncates that tract, so neither mate carries the true repeat length. That
corrupts the evidence indel callers use in repeats, and the microsatellite
jitter models that feed them (e.g. REDUX → SAGE in oncoanalyser).

`--ms_shift` moves the split point by up to N bases to land outside the
tract. Because TLEN and duplicate identity depend only on the read's start
and end, shifting the split changes nothing except the relative mate lengths
— it is free.

Measured on a real quality-trimmed sample: **10.8% of reads had a midpoint
split falling inside a microsatellite, and a shift of ≤6 bp avoided 100% of
them.** On a synthetic repeat-rich set, 6.1% were affected and 89.8% were
avoidable, leaving 0.62% residual. The per-sample counts are reported in
`*.stats.tsv` and in the MultiQC table.

## Effect on downstream analysis

The converted data preserves every base of every retained read, but it is
still derived data. Two properties drive everything below:

1. **There are no discordant pairs, ever.** Every pair is FR with
   TLEN = read length. Orientation and fragment-size anomalies — the evidence
   structural-variant callers are built on — cannot exist by construction.
2. **The two mates are halves of one physical read.** They do not carry
   independent sequencing errors, and a feature longer than a mate (an MNV,
   an indel, a repeat tract, a CDR3) may be divided across the pair.

What is *not* a problem, contrary to what the pre-1.2.0 documentation said:
coverage, duplicate marking, mappability, proper-pair rate, and GC profile
are all intact.

### nf-core/oncoanalyser 3.0

Judged per output of [oncoanalyser 3.0](https://nf-co.re/oncoanalyser/3.0.0/docs/output/).
Process names are as accepted by `--processes_exclude`.

| Output | Effect of converted PE input | Severity | Mitigation |
| ------ | ---------------------------- | -------- | ---------- |
| **ALIGNMENT** (BWA-MEM2) | 100% properly paired, no singletons. MAPQ ≥ 30 in 93.0% of mates vs 92.5% of the SE reads — mate rescue offsets the shorter length. | None | — |
| **REDUX** | Duplicate marking exact (verified). Microsatellite metrics are harmed where a split truncates a repeat tract. | Low | `ms_shift` removed all such splits on real data |
| **BAMTOOLS** | Coverage and flag counts correct. Fragment-length distribution equals the read-length distribution, so there is no library mode — accurate for this data, but it will not look like a normal library. | Low (interpretive) | MultiQC table states this explicitly |
| **SAGE** | SNVs sound: coverage, MAPQ and allele fractions intact, and no double-counting since each base appears in exactly one mate. But every base now sits nearer a mate edge, so SAGE's read-edge penalties bite harder; and MNVs or indels spanning the split are divided across the pair. | Low (SNV), Moderate (MNV/indel) | `gap_size=0` maximises mate length; `ms_shift` protects repeat indels |
| **SAGE_APPEND**, **SAGE_VISUALISER**, **QSEE** | Operate on SAGE output. | inherits SAGE | — |
| **PAVE** | Pure annotation of SAGE calls. | inherits SAGE | — |
| **ESVEE** | No discordant pairs exist; TLEN is capped at the read length so no deletion/duplication fragment signal is possible; split-read evidence is divided across mates. Also crashes on UMI-style read IDs (`ReadIdTrimmer`, unrelated to this pipeline). | **Severe** | None possible — exclude |
| **AMBER** | B-allele frequencies at het SNP loci depend on depth and allele fraction only. | None | — |
| **COBALT** | Per-bin read depth unchanged; with `gap_size=0` no bases are lost, so the GC profile is intact. | None | `gap_size=0` removes the former ~10% base loss |
| **PURPLE** | Depth (COBALT) and BAF (AMBER) inputs are sound, so purity, ploidy and broad copy-number segments are reliable. SV input is not, so breakpoint-precise segment boundaries and focal events degrade. | Moderate | Exclude ESVEE so PURPLE segments on depth rather than fitting invalid breakpoints |
| **LINX** / linxreport | Entirely SV-driven: clusters, fusions and SV-derived drivers. | **Severe** | None — exclude |
| **CHORD** | HRD prediction leans on SV counts, SV size distribution and deletion-with-microhomology features. | **Severe** | None — exclude |
| **TEAL** | Halving reads halves the telomeric repeats per read, so fewer reads clear the repeat-count threshold and telomere length is *systematically underestimated* rather than merely noisy. Telomeric rearrangements are SV-based. | **Severe** | None — a fully telomeric read is one continuous tract, so `ms_shift` cannot help. Exclude |
| **VIRUSINTERPRETER** (+ VIRUSBreakend) | Viral integration sites need discordant and split-read evidence. Viral read *content* partly survives; integration calls do not. | **Severe** (integration) | None — exclude |
| **CUPPA** | Mixes SNV features (sound) with SV and copy-number features (invalid or degraded), so the combined classification is unreliable. | Moderate–Severe | Exclude, or treat the call as unreliable |
| **LILAC** | Mates are contiguous and jointly span the original read, so fragment-level evidence is intact; per-read allele discrimination across distant polymorphic sites is weaker. | Moderate | `gap_size=0` keeps mates contiguous with no lost middle |
| **CIDER** | A split falling inside a CDR3 divides it across the mates; fragment-level assembly recovers part of this. | Moderate | `gap_size=0` |
| **NEO** | Fusion-derived neoepitopes inherit LINX and are invalid. Point-mutation-derived neoepitopes are sound. | Moderate | Exclude if fusion neoepitopes matter |
| **SIGS** | Trinucleotide context is taken from the reference, not the read, and the SNV calls are sound. | Low | — |
| **PEACH** | Germline SNV haplotypes at known pharmacogenomic loci. | None | — |
| **WISP** | Tumour fraction is driven by variant allele fractions and copy number, which are sound. Any fragment-length (fragmentomics) component is destroyed, since TLEN is the read length. | Low–Moderate | — |
| **ISOFOX** | RNA quantification — not applicable to DNA input. | N/A | — |
| **ORANGE** | Aggregates every module, including the invalid ones, and presents them with equal authority in a single PDF. | **Interpretive hazard** | Exclude the invalid modules so they never reach the report |
| **MULTIQC** / Nextflow reports | Run QC and provenance. | None | — |

**Suggested exclusion list:**

```bash
--processes_exclude esvee,linx,chord,teal,virusinterpreter,cuppa,neo,cider
```

Verify that PURPLE still runs with ESVEE excluded — PURPLE consumes ESVEE
output, and depending on the version it may require the input to be present
rather than merely empty. If it refuses, the alternative is to let ESVEE run
and then discard everything downstream of it, rather than to trust its calls.

### Other somatic callers

The same reasoning applies outside oncoanalyser:

- **Reliable**: Mutect2, Strelka2, VarScan2, LoFreq and other SNV/indel
  callers; depth-based CNV callers (CNVkit, Control-FREEC).
- **Do not use**: Manta, Delly, GRIDSS, LUMPY, SvABA or any other PE-aware
  structural-variant caller.
- **Strand-bias filters** (Mutect2 `StrandArtifact`, Strelka2's strand
  filters) are uninformative: every variant appears on both "strands" by
  construction, since R1 is forward and R2 is the reverse of the same
  physical base call. Consider disabling them.
- **Error independence**: a sequencing error appears in exactly one mate, so
  there is no double-counting — but "two reads supporting a variant" across a
  pair still corresponds to one physical observation.

For structural variants, run SE-aware tools on the **original SE FASTQ** in a
separate branch (`pindel -iSE`, GRIDSS `--singleend`, or depth-based callers),
rather than trying to make one converted dataset serve both purposes.

## Credits

Biocentric/nf-se-to-pe is developed and maintained by Biocentric.
Pipeline structure follows [nf-core](https://nf-co.re/) conventions.

## Citations

- SAMtools: Danecek P, et al. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, 10(2): giab008.
- FastQC: Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data.
- MultiQC: Ewels P, et al. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*, 32(19): 3047–3048.
- Nextflow: Di Tommaso P, et al. (2017). Nextflow enables reproducible computational workflows. *Nature Biotechnology*, 35: 316–319.

## License

See [LICENSE](LICENSE).
