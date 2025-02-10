process GENESPACE_RUN {
    tag "$meta.id"
    label 'process_low'

    // conda "${moduleDir}/environment.yml"
    conda "/users/nadjafn/.conda/envs/orthofinder"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hdcf5f25_4' :
        'biocontainers/gffread:0.12.7--hdcf5f25_4' }"

    input:
    tuple val(meta), val(haplotypes), path(fasta), path(gff), path(input_dir)
    path (MCscanX)

    output:
    tuple val(meta), path("*.tsv")  , emit: dir             , optional: true
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args             ?: ''
    def prefix      = task.ext.prefix           ?: "${meta.id}_genespace"
    def output_name = "${prefix}.tsv"
    def ref_haplotype = haplotypes[0]
    """
    Rscript $baseDir/scripts/run_Genespace.R --working_dir $input_dir \
                                             --mcscanx_path $MCscanX  \
                                             --ploidy 1 \
                                             --ref_genome hap1 \
                                             --output $output_name \

    Rscript --version > versions.yml
    """

}

