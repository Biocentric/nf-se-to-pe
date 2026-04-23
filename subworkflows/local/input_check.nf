include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check/main'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header: true, sep: ',' )
        .map { row ->
            def meta = [ id: row.sample, single_end: true ]
            def fq   = file(row.fastq, checkIfExists: true)
            [ meta, fq ]
        }
        .set { reads }

    emit:
    reads    = reads
    versions = SAMPLESHEET_CHECK.out.versions
}
