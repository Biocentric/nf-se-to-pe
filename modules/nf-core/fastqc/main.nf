process FASTQC {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::fastqc=0.12.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0' :
        'biocontainers/fastqc:0.12.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def renamed = reads instanceof List
        ? reads.withIndex().collect { f, i -> "ln -s ${f} ${prefix}_${i + 1}.${f.name.endsWith('.gz') ? 'fastq.gz' : 'fastq'}" }.join('\n    ')
        : "ln -s ${reads} ${prefix}.${reads.name.endsWith('.gz') ? 'fastq.gz' : 'fastq'}"
    """
    ${renamed}
    fastqc ${args} --threads ${task.cpus} ${prefix}*.fastq*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed -e "s/FastQC v//g")
    END_VERSIONS
    """
}
