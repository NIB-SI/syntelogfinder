
process look_at_exact_dups {
    publishDir(
        path: "${params.publish_dir}/Duplicated",
    )
    cache 'lenient'
    conda "/users/nadjafn/.conda/envs/orthofinder"
    input:
        tuple path(blast), val(meta)

    output:
        path("$meta")
              
    script:
    """
    python /scratch/nadjafn/potato-allelic-orthogroups/scripts/cluster_diamon.py ${blast} ${blast}_out.tsv
    mkdir $meta 
    mv ${blast}_out.tsv $meta 
    mv .command.out $meta
    mv *png $meta
    """
}