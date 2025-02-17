process GENESPACE_PARSE {
    tag "$meta.id"
    label 'process_low'

    // conda "${moduleDir}/environment.yml"
    conda "/DKED/scratch/nadjafn/Phasing/ASE-tools-benchmark/conda/expressionMatrix-613c97c23a72e82e6de2bbc3b086d489"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hdcf5f25_4' :
        'biocontainers/gffread:0.12.7--hdcf5f25_4' }"

    input:
    tuple val(meta), val(haplotypes), path(pangenes), path(gff)
  

    output:
    tuple val(meta), val(haplotypes), path("*.tsv")  , emit: pangenes
    tuple val(meta), path("*.gff")  , emit: gff
    tuple val(meta), path("*exploded*"), emit: syntelogs
    tuple val(meta), path("*.png")  , emit:plots
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args             ?: ''
    def prefix      = task.ext.prefix           ?: "${meta.id}_genespace"

    """
    python $baseDir/scripts/parse_pangenes.py --pangenes $pangenes \
                                              --gff $gff \
                                              --output $prefix                         
    python --version > versions.yml
    """

}

