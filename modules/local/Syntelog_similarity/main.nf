#!/usr/bin/env nextflow

process SYNTELOG_SIMILARITY {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/pip_argparse_gffutils_matplot_pruned:3b26a631aedc4256' :
        'community.wave.seqera.io/library/pip_argparse_gffutils_matplot_pruned:3b26a631aedc4256' }"

    input:
        tuple val(meta), val(haplotypes), path(pangenes)
        tuple val(meta), path(blast_file)

    output:
        tuple val(meta), val(haplotypes), path("${prefix}*_analysis.tsv"), emit: blast_pangenes
        tuple val(meta), path("${prefix}*.svg"), emit: plots
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        prefix = task.ext.prefix ?: "${meta.id}_blast"


    """
    CDS_similarity_BLAST.py \\
        --blast ${blast_file} \\
        --syntenic_genes ${pangenes} \\
        --output ${prefix} \\
        --ploidy ${params.ploidy} \\
        $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        CDS_similarity_BLAST: \$(python $projectDir/bin/CDS_similarity_BLAST.py --version 2>&1 | sed 's/CDS_similarity_BLAST.py //g')
    END_VERSIONS
    """
    stub:
    prefix = task.ext.prefix ?: "${meta.id}_blast"
    """
    touch ${prefix}.tsv
    touch ${prefix}_plot1.png
    touch ${prefix}_plot2.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        CDS_similarity_BLAST: \$(python $projectDir/scripts/CDS_similarity_BLAST.py --version 2>&1 | sed 's/CDS_similarity_BLAST.py //g')
    END_VERSIONS
    """
}
