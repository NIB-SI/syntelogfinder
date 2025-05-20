include { GFF_TO_PROMOTER } from '../../modules/local/gff_to_promotor'


workflow PROMOTOR_COMPARISON {
    take:
        gff_file
        fasta_ch
        promotor_length

    main:
        GFF_TO_PROMOTER(gff_file, fasta_ch, promotor_length)

    emit:
        genespace_parse.pangenes
}