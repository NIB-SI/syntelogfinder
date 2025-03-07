process GENESPACE_PARSE {
    tag "$meta.id"
    label 'process_low'

    // conda "${moduleDir}/environment.yml"
    conda "/DKED/scratch/nadjafn/Phasing/ASE-tools-benchmark/conda/expressionMatrix-613c97c23a72e82e6de2bbc3b086d489"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pandas_python_pip_argpar_pruned:72e2ed5052546765':
        'biocontainers/gffread:0.12.7--hdcf5f25_4' }"

    input:
    tuple val(meta), val(haplotypes), path(pangenes), path(gff)
    val(cdsLengthExtension)
  

    output:
    tuple val(meta), val(haplotypes), path("*.tsv")  , emit: pangenes
    tuple val(meta), path("*.png")  , emit:plots
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args             ?: ''
    def prefix      = task.ext.prefix           ?: "${meta.id}_genespace"

    """
    python $baseDir/scripts/parse_genespace_pangenes.py --pangenes $pangenes \
                                                        --gff $gff \
                                                        --extend_CDS $cdsLengthExtension \
                                                        --output $prefix \
                                                        -s 1hap1_1hap2_1hap3_1hap4_synteny                 
    python --version > versions.yml

}

