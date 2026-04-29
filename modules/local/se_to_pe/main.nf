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
    path "versions.yml",                                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def gap        = params.gap_size     ?: 10
    def min_mate   = params.min_mate_len ?: 30
    def jitter_max = params.jitter_max   ?: 0
    def seed       = params.seed         ?: 42
    """
    split_se_to_pe.py \\
        --in ${reads} \\
        --r1 ${prefix}_1.fastq.gz \\
        --r2 ${prefix}_2.fastq.gz \\
        --gap-size ${gap} \\
        --jitter-max ${jitter_max} \\
        --seed ${seed} \\
        --min-mate-len ${min_mate} \\
        --stats ${prefix}.stats.tsv

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
    echo -e "reads_in\\treads_out\\treads_dropped\\tbp_in\\tbp_out\\tbp_retained_pct\\n0\\t0\\t0\\t0\\t0\\t0.00" > ${prefix}.stats.tsv
    touch versions.yml
    """
}
