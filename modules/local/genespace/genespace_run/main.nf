process GENESPACE_RUN {
    tag "$meta.id"
    label 'process_medium'

    cache 'lenient'

    conda "${moduleDir}/environment.yml"
    //conda "/users/nadjafn/.conda/envs/genespace-env"
    // https://github.com/HuffordLab-Containers/genespace_docker/blob/main/Dockerfile
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://arnstrm2/genespace:1.1.4' :
        'arnstrm2/genespace:1.1.4' }"
    // Problem optparse and library(dplyr) not in conatiner


    input:
    tuple val(meta), val(haplotypes), path(fasta), path(gff), path(input_dir)

    output:
    tuple val(meta), val(haplotypes), path("*.tsv")  , emit: pangenes      , optional: true
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args             ?: ''
    def prefix      = task.ext.prefix           ?: "${meta.id}_genespace"
    def output_name = "${prefix}.tsv"
    def ref_haplotype = haplotypes[1]

    """
    run_Genespace.R --working_dir $input_dir \
                                             --mcscanx_path $params.mcscanx_path \
                                             --ploidy 1 \
                                             --ref_genome $ref_haplotype \
                                             --threads $task.cpus \
                                             --output $output_name

    Rscript --version > versions.yml
    """

}
