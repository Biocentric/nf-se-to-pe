# ![nf-se-to-pe](docs/images/nf-se-to-pe_logo_light.png)

[![GitHub Actions CI Status](https://github.com/Biocentric/nf-se-to-pe/actions/workflows/ci.yml/badge.svg)](https://github.com/Biocentric/nf-se-to-pe/actions/workflows/ci.yml)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed.svg?labelColor=000000)](https://www.docker.com/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049.svg?labelColor=000000)](https://docs.conda.io/en/latest/)

## Introduction

**Biocentric/nf-se-to-pe** is a bioinformatics pipeline that converts
single-end (SE) FASTQ data into simulated paired-end (PE) FASTQ data so
that downstream tools requiring PE input (e.g. somatic variant callers
that expect `_R1`/`_R2` files) can accept it.

The conversion is **sequence-lossless-by-construction**: every output base
comes from the original SE read, no synthetic sequence or error noise is
introduced. Each SE read is split near the midpoint; R2 is the
reverse-complement of the second half. A configurable gap between R1 and
R2 simulates an insert. An optional small per-read random jitter on the R2
start position breaks the otherwise-degenerate fragment-endpoint
distribution that would cause downstream tools (notably
`MarkDuplicates`) to flag all pairs as duplicates.

## Pipeline summary

1. Samplesheet validation ([`check_samplesheet.py`](bin/check_samplesheet.py))
2. Raw FASTQ QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
3. SE → PE conversion ([`split_se_to_pe.py`](bin/split_se_to_pe.py))
4. Converted FASTQ QC (`FastQC`)
5. Aggregate QC report ([`MultiQC`](https://multiqc.info/))

## Quick start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=23.04.0`).

2. Install one of the container/conda runtimes: [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), or [`Conda`](https://conda.io/miniconda.html).

3. Test the pipeline on a minimal dataset:

    ```bash
    nextflow run Biocentric/nf-se-to-pe -profile test,singularity --outdir results_test
    ```

4. Run on your own data:

    ```bash
    nextflow run Biocentric/nf-se-to-pe \
        -profile singularity \
        --input samplesheet.csv \
        --outdir results \
        --gap_size 10 \
        --jitter_max 3
    ```

## Samplesheet input

You need to create a CSV samplesheet with your input single-end FASTQ
files. It should look like this:

```csv
sample,fastq
tumor_01,/path/to/tumor_01.fastq.gz
normal_01,/path/to/normal_01.fastq.gz
```

| Column   | Description                                                                                                              |
| -------- | ------------------------------------------------------------------------------------------------------------------------ |
| `sample` | Custom sample name. Must be unique. Spaces are not allowed.                                                              |
| `fastq`  | Full path or URL to a single-end FASTQ file. File must end in `.fastq`, `.fq`, `.fastq.gz`, or `.fq.gz`.                 |

## Parameters

| Parameter          | Description                                                                                              | Default |
| ------------------ | -------------------------------------------------------------------------------------------------------- | ------- |
| `--input`          | Path to the input samplesheet CSV. **Required.**                                                         | `null`  |
| `--outdir`         | The output directory where the results will be saved.                                                    | `./results` |
| `--gap_size`       | Number of bases dropped between R1 and R2 to simulate insert gap.                                        | `10`    |
| `--jitter_max`     | Maximum extra bases dropped at R2 start per read, drawn uniformly from `[0, jitter_max]`. Set `0` to disable. | `3`     |
| `--seed`           | RNG seed for jitter (reproducibility).                                                                   | `42`    |
| `--min_mate_len`   | Minimum length required per output mate. Shorter input reads are dropped.                                | `30`    |
| `--singularity_cache_dir` | Override the Singularity image cache directory.                                                   | `~/.singularity-cache/se-to-pe` |

## How the conversion works

For an input SE read of length `L`, the conversion is **adaptive** — gap
and jitter shrink for short reads so that data with variable read
lengths (e.g. quality-trimmed reads) is preserved as much as possible:

- `max_extra = L − 2 × min_mate_len`
- Reads shorter than `2 × min_mate_len` are dropped
- `effective_gap = min(gap_size, max_extra)`
- `j ~ Uniform(0, min(jitter_max, max_extra − effective_gap))`
- `drop = effective_gap + j`, `usable = L − drop`
- `r1_len = ⌈usable / 2⌉`, `r2_len = usable − r1_len`
- **R1** = `seq[0 : r1_len]`, original quality
- **R2** = reverse-complement of `seq[r1_len + drop : r1_len + drop + r2_len]`, with quality reversed (not complemented)

For long reads the configured `gap_size` and `jitter_max` are used in
full. For reads that can only fit two minimum-length mates, gap and
jitter shrink to 0 and the read is split exactly in half.

Read headers are rewritten to pair cleanly: Illumina-style headers with
`1:N:0:…` become `1:N:0:…` / `2:N:0:…`; `/1` suffixes become `/1` / `/2`;
bare headers get `/1` / `/2` appended.

## Downstream caveats for somatic variant calling

The generated PE data preserves the original bases exactly but is still
semi-synthetic. Important things to know:

1. **MarkDuplicates**: With `jitter_max ≥ 3`, R2 reference-end position
   varies enough that genuine PCR duplicate clusters are *not* fully
   re-grouped. You can run MarkDuplicates, but it will under-call
   duplicates. Increase `--jitter_max` (e.g. to 7) to suppress residual
   dup-grouping further — at the cost of a few more dropped bases per
   read.
2. **Strand bias filters** (e.g. Mutect2 `StrandArtifact`,
   `STRAND_BIAS`, Strelka2's strand filters): every variant appears on
   both "strands" by construction (R1 forward, R2 reverse of the same
   physical base-call). These filters become uninformative — consider
   disabling them.
3. **Insert-size distribution**: narrowly distributed around `gap_size`.
   Do not tune insert-size-dependent parameters based on what you see.
4. **Error independence**: a sequencing error on the original read
   appears on both R1 and R2 at mirrored positions. "Two reads supporting
   a variant" in the converted data corresponds to one physical
   observation.

## Credits

Biocentric/nf-se-to-pe is developed and maintained by Biocentric.

Pipeline structure follows [nf-core](https://nf-co.re/) conventions.

## Citations

If you use this pipeline, please cite the tools it wraps:

- FastQC: Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data.
- MultiQC: Ewels P, et al. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*, 32(19): 3047–3048.
- Nextflow: Di Tommaso P, et al. (2017). Nextflow enables reproducible computational workflows. *Nature Biotechnology*, 35: 316–319.

## License

See [LICENSE](LICENSE).
