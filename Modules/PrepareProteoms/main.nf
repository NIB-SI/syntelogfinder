
process download_rename {
    cache 'lenient'
    input:
        val(cultivar)

    output:
        path("GFF_Split")
        path("Haplotype_Split")
        path('cDNA_Split')
              
    script:
    """
    bash $baseDir/scripts/download_and_rename.sh
    """
}


process merge_gff_to_bed {
    cache 'lenient'
    conda "/users/nadjafn/.conda/envs/orthofinder"
    input:
        path(gff_dir)

    output:
        path("cultivars.bed")

    script:
    """
    cat $gff_dir/* | awk -v OFS='\\t' '\$3 == "mRNA" { split(\$9, a, ";"); split(a[1], b, "="); print \$1, b[2], \$4, \$5 }' > cultivars.bed
    """
}


process get_CDS {
    cache 'lenient'
    conda "/users/nadjafn/.conda/envs/gffread"
    input:
        tuple val(culitvar),
              path(fasta),
              path(gff)

    output:
        tuple val(culitvar),
              path("${culitvar}.pep")
              
    script:
    """
    gffread -y ${culitvar}.pep -g ${fasta} ${gff}
    """
}



process split_into_haplotypes {
    cache 'lenient'
    
    input:
        tuple val(culitvar),
              path(fasta),
              path(gff)

    output:
        tuple val(culitvar),
              path("${culitvar}.pep")
              
    script:
    """
    gffread -y ${culitvar}.pep -g ${fasta} ${gff}
    """
}

