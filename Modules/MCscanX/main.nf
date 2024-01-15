process run_MCScanX {
    cache 'lenient'
    conda "/users/nadjafn/.conda/envs/orthofinder"
    input:
        path(merged_bed)
        path(blast_output)

    output:
        path("scan_input/cultivar*")
              
    script:
    """
    mkdir scan_input 
    cp $merged_bed scan_input/cultivar.gff
    cp $blast_output scan_input/cultivar.blast
    /scratch/nadjafn/MCScanX/MCScanX scan_input/cultivar 
    """
}
