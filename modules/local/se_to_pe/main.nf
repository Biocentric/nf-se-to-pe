process SE_TO_PE {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.11"
    container 'python:3.11'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_1.fastq.gz"), path("*_2.fastq.gz"), emit: reads
    tuple val(meta), path("*.stats.tsv"),                        emit: stats
    tuple val(meta), path("*_mqc.tsv"),                          emit: mqc
    path "versions.yml",                                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Use explicit null checks, not Groovy's elvis operator: 0 is falsy, so
    // `params.gap_size ?: 10` would silently turn a deliberate 0 into 10.
    def gap        = params.gap_size     != null ? params.gap_size     : 0
    def jitter_max = params.jitter_max   != null ? params.jitter_max   : 0
    def ms_shift   = params.ms_shift     != null ? params.ms_shift     : 6
    def min_mate   = params.min_mate_len != null ? params.min_mate_len : 30
    def seed       = params.seed         != null ? params.seed         : 42
    """
    split_se_to_pe.py \\
        --in ${reads} \\
        --r1 ${prefix}_1.fastq.gz \\
        --r2 ${prefix}_2.fastq.gz \\
        --gap-size ${gap} \\
        --jitter-max ${jitter_max} \\
        --ms-shift ${ms_shift} \\
        --seed ${seed} \\
        --min-mate-len ${min_mate} \\
        --sample ${prefix} \\
        --stats ${prefix}.stats.tsv \\
        --mqc ${prefix}.se_to_pe_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_1.fastq.gz
    echo "" | gzip > ${prefix}_2.fastq.gz
    printf 'reads_in\\treads_out\\treads_dropped\\tbp_in\\tbp_out\\tbp_retained_pct\\tms_splits_shifted\\tms_splits_residual\\n0\\t0\\t0\\t0\\t0\\t0.00\\t0\\t0\\n' > ${prefix}.stats.tsv
    printf '# id: se_to_pe\\nSample\\tReads in\\n${prefix}\\t0\\n' > ${prefix}.se_to_pe_mqc.tsv
    touch versions.yml
    """
}
