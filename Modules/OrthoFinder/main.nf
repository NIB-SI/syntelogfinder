process run_OrthoFinder {
    publishDir(
        path: "${params.publish_dir}/OrthoFinder",
    )
    cache 'lenient'
    conda "/users/nadjafn/.conda/envs/orthofinder"
    input:
        path(prot_dir)
        val(prefix)

    output:
        path("Orthofinder_output")
    script:
    """
    mkdir filtered_prot_dir 
    cp ${prot_dir}/${prefix}* filtered_prot_dir 
    orthofinder -t 60 -f filtered_prot_dir -o Orthofinder_output
    """
}