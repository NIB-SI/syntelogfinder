include { GFF_TO_PROMOTER } from '../../modules/local/gff_to_promotor'


workflow PROMOTOR_COMPARISON {
    take:
        gff_file
        fasta_ch
        promotor_length

    main:
        promotor_gff = GFF_TO_PROMOTER(gff_file, fasta_ch, promotor_length)

    // emit:
    //     promotor_gff
}