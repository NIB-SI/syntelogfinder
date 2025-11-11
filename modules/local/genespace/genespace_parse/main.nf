process GENESPACE_PARSE {
    tag "$meta.id"
    label 'process_low'

    cache 'lenient'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pip_argparse_gffutils_matplotlib_pruned:acaa38c7b16e2893' :
        'community.wave.seqera.io/library/pip_argparse_gffutils_matplotlib_pruned:3c48a85c80b84759' }"


    input:
    tuple val(meta), val(haplotypes), path(pangenes), path(gff)

    output:
    tuple val(meta), val(haplotypes), path("${prefix}_categories.tsv"), emit: pangenes
    tuple val(meta), path("${prefix}*.svg"), emit: plots
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_genespace"
    def haplotypes_arg = haplotypes.sort().collect { "1${it}" }.join('_') + '_s'

    """
    parse_genespace_pangenes.py \\
        --pangenes $pangenes \\
        --gff $gff \\
        --output ${prefix} \\
        -s $haplotypes_arg \\
        $args
    echo "echo!!?!?##!"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        parse_genespace_pangenes: \$(python $projectDir/bin/parse_genespace_pangenes.py --version 2>&1 | sed 's/parse_genespace_pangenes.py //g')
    END_VERSIONS
    """
    stub:
    prefix = task.ext.prefix ?: "${meta.id}_genespace"
    """
    touch ${prefix}.tsv
    touch ${prefix}_plot1.png
    touch ${prefix}_plot2.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        parse_genespace_pangenes: \$(python $projectDir/scripts/parse_genespace_pangenes.py --version 2>&1 | sed 's/parse_genespace_pangenes.py //g')
    END_VERSIONS
    """
}
