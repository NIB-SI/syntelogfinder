// Run diamond to align protein sequences to all

process run_blastn {
    cache 'lenient'
    conda "/users/nadjafn/.conda/envs/blast"
    input:
        path(cDNA_dir)
        val(prefix)
            
    output:
        tuple path("nucleotide_cultivar.blast"), val('nuc_blast')
              
    script:
    """
    cat ${cDNA_dir}/${prefix}* > all_cdna.fa
    makeblastdb -dbtype nucl -out all_cdna_db.fa -in all_cdna.fa 
    blastn -db all_cdna_db.fa  -query all_cdna.fa  -out nucleotide_cultivar.blast -perc_identity 100  -qcov_hsp_perc 100 -outfmt 6
    """
}