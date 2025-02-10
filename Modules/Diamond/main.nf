// Run diamond to align protein sequences to all

process DIAMOND {
    cache 'lenient'
    conda "/users/nadjafn/.conda/envs/orthofinder"
    input:
        tuple val(meta), path(pep)
            
    output:
        tuple val(meta), path("${meta.id}.blast")
              
    script:
    """
    diamond makedb --in $pep -d ${meta.id}_db
    diamond blastp -p 30 -d ${meta.id}_db -q $pep -o ${meta.id}.blast --max-target-seqs 20 --evalue 0.000001 --no-self-hits
    """
}


// Run diamond to align protein sequences to all
process filter_diamond {
    publishDir(
        path: "${params.publish_dir}/Duplicated",
    )
    cache 'lenient'
    conda "/users/nadjafn/.conda/envs/orthofinder"
    input:
        path(prot_dir)
        val(prefix)
            
    output:
        tuple path("cultivar.blast"), val('prot_diamond')
              
    script:
    """
    cat ${prot_dir}/${prefix}* > all_prot.fa
    diamond makedb --in all_prot.fa --db all_prot_db.fa
    diamond blastp -p 30 -d all_prot_db.fa.dmnd -q all_prot.fa  -o cultivar.blast --max-target-seqs 20 --evalue 0.000001 --id 100 --query-cover 100
    python /scratch/nadjafn/potato-allelic-orthogroups/scripts/cluster_diamon.py cultivar.blast  out.tsv
    """
}



// get uniuq
// awk -F"\t" '$1 != $2' cultivar.cluster |wc -l