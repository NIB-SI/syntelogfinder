include { GENESPACE_INPUT_PREPERATION } from '../../modules/local/genespace/genespace_input_preperation'
include { GENESPACE_RUN } from '../../modules/local/genespace/genespace_run'
include { GENESPACE_PARSE } from '../../modules/local/genespace/genespace_parse'

workflow GENESPACE_ANALYSIS {
    take:
        gffread_output
        mcscanx_path
        agat_output

    main:
        genespace_input = GENESPACE_INPUT_PREPERATION(gffread_output)
        genespace_run = GENESPACE_RUN(genespace_input.dir, mcscanx_path)
        genespace_parse = GENESPACE_PARSE(genespace_run.pangenes.join(agat_output))

    emit:
        genespace_parse.pangenes
}