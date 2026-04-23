include { INPUT_CHECK       } from '../subworkflows/local/input_check'
include { SE_TO_PE          } from '../modules/local/se_to_pe/main'
include { FASTQC as FASTQC_IN  } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_OUT } from '../modules/nf-core/fastqc/main'
include { MULTIQC           } from '../modules/nf-core/multiqc/main'

workflow SE_TO_PE_PIPELINE {

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    INPUT_CHECK ( file(params.input, checkIfExists: true) )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    FASTQC_IN ( INPUT_CHECK.out.reads )
    ch_versions = ch_versions.mix(FASTQC_IN.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_IN.out.zip.collect{ it[1] })

    SE_TO_PE ( INPUT_CHECK.out.reads )
    ch_versions = ch_versions.mix(SE_TO_PE.out.versions.first())

    ch_pe = SE_TO_PE.out.reads.map { meta, r1, r2 ->
        [ meta + [single_end: false], [ r1, r2 ] ]
    }
    FASTQC_OUT ( ch_pe )
    ch_versions = ch_versions.mix(FASTQC_OUT.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_OUT.out.zip.collect{ it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(SE_TO_PE.out.stats.collect{ it[1] })

    ch_multiqc_config = params.multiqc_config
        ? Channel.fromPath(params.multiqc_config, checkIfExists: true)
        : Channel.empty()
    ch_multiqc_logo   = Channel.empty()
    ch_multiqc_extra  = Channel.empty()

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_extra.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    reads    = SE_TO_PE.out.reads
    stats    = SE_TO_PE.out.stats
    multiqc  = MULTIQC.out.report
    versions = ch_versions
}
