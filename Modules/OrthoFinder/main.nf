process run_OrthoFinder {
    publishDir(
        path: "${params.publish_dir}/OrthoFinder",
    )
    cache 'lenient'
    conda "/users/nadjafn/.conda/envs/orthofinder"
    input:
        path(prot_dir)

    output:
        path("Orthofinder_output")
    script:
    """
    orthofinder -t 60 -f ${prot_dir}  -o Orthofinder_output
    """
}