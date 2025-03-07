process GENESPACE_PARSE {
    tag "$meta.id"
    label 'process_low'

    conda "/DKED/scratch/nadjafn/Phasing/ASE-tools-benchmark/conda/expressionMatrix-613c97c23a72e82e6de2bbc3b086d489"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pandas_python_pip_argpar_pruned:72e2ed5052546765':
        'biocontainers/gffread:0.12.7--hdcf5f25_4' }"

    input:
    tuple val(meta), val(haplotypes), path(pangenes), path(gff)

    output:
    tuple val(meta), val(haplotypes), path("${prefix}.tsv"), emit: pangenes
    tuple val(meta), path("${prefix}*.png"), emit: plots
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_genespace"
    def haplotypes_arg = haplotypes.join('_') + '_synteny'

    """
    python $projectDir/bin/parse_genespace_pangenes.py \\
        --pangenes $pangenes \\
        --gff $gff \\
        --output ${prefix} \\
        -s ${haplotypes_arg} \\
        $args

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
        parse_genespace_pangenes: \$(python $projectDir/bin/parse_genespace_pangenes.py --version 2>&1 | sed 's/parse_genespace_pangenes.py //g')
    END_VERSIONS
    """
}