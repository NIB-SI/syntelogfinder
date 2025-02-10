process run_MCScanX {
    cache 'lenient'
    conda "/users/nadjafn/.conda/envs/orthofinder"
    input:
        tuple val(meta1), path(bed)
        tuple val(meta2), path(blast_output)

    output:
        path("scan_input/cultivar*")
              
    script:
    """
    bash /scratch/nadjafn/potato-allelic-orthogroups/scripts/renameBED.sh  $bed ${bed}_renamed

    
    mkdir scan_input 
    cp ${bed}_renamed scan_input/cultivar.gff
    cp $blast_output scan_input/cultivar.blast
    /scratch/nadjafn/MCScanX/MCScanX scan_input/cultivar 
    """
}