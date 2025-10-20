#!/usr/bin/env nextflow

process SYNTELOG_SIMILARITY {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pandas_python_pip_argpar_pruned:72e2ed5052546765':
        'biocontainers/gffread:0.12.7--hdcf5f25_4' }"
    input:
    tuple val(meta), val(haplotypes), path(pangenes)
    tuple val(meta), path(blast_file)

    output:
    tuple val(meta), val(haplotypes), path("${prefix}*_analysis.tsv"), emit: blast_pangenes
    tuple val(meta), path("${prefix}*.png"), emit: plots
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_blast"


    """
    python $projectDir/scripts/CDS_similarity_BLAST.py \\
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
