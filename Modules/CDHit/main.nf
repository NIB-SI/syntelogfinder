
// Run cs-hit-est to find exact gene dublications

process run_cd_hit {
    cache 'lenient'
    conda "/users/nadjafn/.conda/envs/orthofinder"
    input:
        path(prot_dir)
        val(prefix)
            
    output:
        path("cd-hit.out.clstr")
              
    script:
    """
    cat ${prot_dir}/${prefix}* > all_prot.fa
    cd-hit-est -i all_prot.fa -o cd-hit.out -sc 1 -d 27 -aS 1 -s 1
    python $baseDir/scripts/parse_CDHit.py cd-hit.out.clstr
    """
}



// only same length  -A 1

// -uS percentage unmatched of shorter seuqence
// cd-hit-est -i all_prot.fa -o cd-hit.out -d 27 -sc 1 -aS 1 -G 0 -g 1 -uS 10 -s 0.8
