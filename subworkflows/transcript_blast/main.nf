include { GFFREAD as GFFREAD_TRANSCRIPT }  from '../../modules/nf-core/gffread'
include { BLAST_MAKEBLASTDB             } from '../../modules/nf-core/blast/makeblastdb'
include { BLAST_BLASTN                  } from '../../modules/nf-core/blast/blastn'


workflow TRANSCRIPT_BLAST {
    take:
        gff             // GFF file from genespace parse
        fasta          // Input FASTA file

    main:
        // Extract transcript sequences
        GFFREAD_TRANSCRIPT (
            gff,
            fasta
        )

        // Build BLAST database
        BLAST_MAKEBLASTDB(
            GFFREAD_TRANSCRIPT.out.gffread_fasta
        )

        // Run BLAST analysis
        BLAST_BLASTN(
            GFFREAD_TRANSCRIPT.out.gffread_fasta,
            BLAST_MAKEBLASTDB.out.db,
            [],                                     // path taxidlist (empty)
            [],                                     // val taxids (empty)
            []
        )

    emit:
        cds_fasta = GFFREAD_TRANSCRIPT.out.gffread_fasta
        results   = BLAST_BLASTN.out.txt
}
