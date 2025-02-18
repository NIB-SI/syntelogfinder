workflow CDS_BLAST {
    take:
        gff             // GFF file from genespace parse
        fasta          // Input FASTA file
    
    main:
        // Extract CDS sequences
        gffread_out = GFFREAD_CDS(
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
        results   = blast_results
}