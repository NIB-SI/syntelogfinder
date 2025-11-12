process GENESPACE_RUN {
    tag "$meta.id"
    label 'process_medium'
    cache 'lenient'
    conda "${moduleDir}/environment.yml"
//conda "/users/nadjafn/.conda/envs/genespace-env"
// https://github.com/HuffordLab-Containers/genespace_docker/blob/main/Dockerfile
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://aewebb/genespace:20250801' :
        'aewebb/genespace:20250801'}"
// Problem optparse and library(dplyr) not in container

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

    // Set MCScanX path based on container engine
    def mcscanx_path = workflow.containerEngine == 'singularity' ?
        '/opt/conda/envs/genespace/bin' :
        params.mcscanx_path

    """
    run_Genespace.R --working_dir $input_dir \\
                    --mcscanx_path $mcscanx_path \\
                    --ploidy 1 \\
                    --ref_genome $ref_haplotype \\
                    --threads $task.cpus \\
                    --output $output_name

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(Rscript --version 2>&1 | sed -n 's/R scripting front-end version //p')
        GENESPACE: \$(Rscript -e "cat(as.character(packageVersion('GENESPACE')))" 2>/dev/null)
    END_VERSIONS
    """
}
