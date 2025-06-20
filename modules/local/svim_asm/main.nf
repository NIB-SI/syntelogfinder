process SVIM_ASM {
    tag "${meta.id.id}_${meta.synt_id}_${meta.direction}"
    label 'process_medium'

    conda "/users/nadjafn/.conda/envs/syri"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/svim-asm:1.0.3--pyhdfd78af_0' :
        'quay.io/biocontainers/svim-asm:1.0.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(query), path(reference)

    output:
    tuple val(meta), path("${prefix}_variants.vcf"), emit: vcf
    tuple val(meta), path("${prefix}"), emit: svim_results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id.id}_${meta.original_query}_${meta.original_reference}"
    

    """
    # Run SVIM-ASM
    svim-asm haploid \\
        $args \\
        ${prefix} \\
        ${bam} \\
        ${reference} \\
        --query_names \\
        --min_sv_size 20

    mv ${prefix}/variants.vcf ${prefix}_variants.vcf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svim-asm: \$(svim-asm --version | sed -e "s/SVIM-ASM v//g")
    END_VERSIONS
    """
}