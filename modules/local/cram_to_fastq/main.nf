process CRAM_TO_FASTQ {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::samtools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(aligned)
    path fasta

    output:
    tuple val(meta), path("*.se.fastq.gz"), emit: reads
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def ref_arg   = fasta ? "--reference ${fasta}" : ''
    def ref_path  = params.ref_path ?: 'https://www.ebi.ac.uk/ena/cram/md5/%s'
    def cache_env = params.ref_cache
        ? "mkdir -p '${params.ref_cache}'\n    export REF_CACHE='${params.ref_cache}/%2s/%2s/%s'"
        : ''
    """
    # CRAM is reference-compressed: htslib needs the reference used at
    # alignment time to decode bases. It resolves this by the M5 checksum in
    # each @SQ header line, checking REF_CACHE first, then REF_PATH. The
    # default REF_PATH points at the EBI/ENA reference registry, so a CRAM
    # aligned to a standard genome decodes without a local FASTA. Set here
    # explicitly rather than relying on an htslib-version-dependent default.
    export REF_PATH='${ref_path}'
    ${cache_env}

    # -F 0x900 drops secondary (0x100) and supplementary (0x800) alignments:
    # those re-present sequence already carried by the primary record, so
    # keeping them would duplicate bases in the output FASTQ.
    #
    # Reads aligned to the reverse strand are stored reverse-complemented in
    # the CRAM; samtools fastq restores original read orientation (and
    # un-reverses quality) automatically for those records.
    #
    # Streams are split by pairing flag so a wrongly-supplied paired-end file
    # is caught rather than silently mangled: genuine SE reads (0x1 unset)
    # land in -0, and -1/-2/-s stay empty.
    samtools fastq \\
        -F 0x900 \\
        ${ref_arg} \\
        --threads ${task.cpus} \\
        -0 ${prefix}.se.fastq.gz \\
        -1 ${prefix}.paired_1.fastq.gz \\
        -2 ${prefix}.paired_2.fastq.gz \\
        -s ${prefix}.singleton.fastq.gz \\
        ${aligned}

    # gzip -cd over several files concatenates them; head -c 1 short-circuits
    # so this stays O(1) regardless of file size. The subshell disables
    # pipefail locally so head's SIGPIPE does not abort the script.
    n_paired=\$( set +o pipefail; gzip -cd \\
        ${prefix}.paired_1.fastq.gz \\
        ${prefix}.paired_2.fastq.gz \\
        ${prefix}.singleton.fastq.gz 2>/dev/null | head -c 1 | wc -c )

    if [ "\$n_paired" -gt 0 ]; then
        echo "ERROR: ${aligned} contains paired-end reads (records with the 0x1 flag set)." >&2
        echo "       nf-se-to-pe converts SINGLE-end data to paired-end; it cannot take" >&2
        echo "       paired-end input. If this file is already paired-end, feed it to your" >&2
        echo "       downstream caller directly instead of running this pipeline." >&2
        exit 1
    fi

    n_se=\$( set +o pipefail; gzip -cd ${prefix}.se.fastq.gz 2>/dev/null | head -c 1 | wc -c )
    if [ "\$n_se" -eq 0 ]; then
        echo "ERROR: no single-end reads extracted from ${aligned}." >&2
        echo "       The file decoded successfully but contained no primary unpaired records." >&2
        exit 1
    fi

    rm -f ${prefix}.paired_1.fastq.gz ${prefix}.paired_2.fastq.gz ${prefix}.singleton.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | sed -n '1s/^samtools //p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.se.fastq.gz
    touch versions.yml
    """
}
