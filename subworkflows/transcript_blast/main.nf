include { GFFREAD as GFFREAD_TRANSCRIPT }  from '../../modules/nf-core/gffread'
include { BLAST_MAKEBLASTDB             } from '../../modules/nf-core/blast/makeblastdb'
include { BLAST_BLASTN                  } from '../../modules/nf-core/blast/blastn'


workflow TRANSCRIPT_BLAST {
    take:
        gff             // GFF file from genespace parse
        fasta          // Input FASTA file

    main:
        // Extract transcript sequences
        gffread_out = GFFREAD_TRANSCRIPT (
            gff,
            fasta
        )

        // Build BLAST database
        blast_db = BLAST_MAKEBLASTDB(
            gffread_out.gffread_fasta
        )

        // Run BLAST analysis
        blast_results = BLAST_BLASTN(
            gffread_out.gffread_fasta,
            blast_db.db
        )

    emit:
        cds_fasta = gffread_out.gffread_fasta
        blast_db  = blast_db.db
        results   = blast_results.txt
}
