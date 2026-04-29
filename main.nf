#!/usr/bin/env nextflow

include { SE_TO_PE_PIPELINE } from './workflows/se_to_pe'

workflow {
    SE_TO_PE_PIPELINE ()
}
