#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SE_TO_PE_PIPELINE } from './workflows/se_to_pe'

workflow {
    SE_TO_PE_PIPELINE ()
}

workflow.onComplete {
    log.info ( workflow.success
        ? "\n[se-to-pe] Pipeline completed successfully. Output: ${params.outdir}"
        : "\n[se-to-pe] Pipeline failed. See ${workflow.workDir} for logs." )
}
