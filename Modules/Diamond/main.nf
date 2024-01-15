// Run diamond to align protein sequences to all

process run_diamond_all {
    cache 'lenient'
    conda "/users/nadjafn/.conda/envs/orthofinder"
    input:
        path(prot_dir)
            
    output:
        path("cultivar.blast")
              
    script:
    """
    cat ${prot_dir}/* > all_prot.fa
    diamond makedb --in all_prot.fa --db all_prot_db.fa
    diamond blastp -p 30 -d all_prot_db.fa.dmnd -q all_prot.fa  -o cultivar.blast   --max-target-seqs 20 --evalue 0.000001
    """
}