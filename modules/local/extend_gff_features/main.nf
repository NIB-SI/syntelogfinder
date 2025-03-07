// modules/local/extend_gff_features.nf

process EXTEND_GFF_FEATURES {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pip_pyranges:b8b5027da61f0906' :
        'quay.io/biocontainers/python:3.8' }"

    input:
    tuple val(meta), path(gff)
    path fasta

    output:
    tuple val(meta), path("*.extended.gff"), emit: gff
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    extend_gff_features.py \\
        --gff $gff \\
        --fasta $fasta \\
        --feature-type CDS \\
        --extension 150 \\
        --output ${prefix}.extended.gff \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        extend_gff_features: \$(extend_gff_features.py --version 2>&1 | sed 's/extend_gff_features.py //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.extended.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        extend_gff_features: \$(extend_gff_features.py --version 2>&1 | sed 's/extend_gff_features.py //g')
    END_VERSIONS
    """
}