// modules/gff_to_promoter.nf

process GFF_TO_PROMOTER {
    tag "$meta.id"
    label 'process_low'

    conda "syri"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/pip_numpy_pyfaidx_pyranges:d084832c91636c20' :
        'quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    tuple val(meta), path(gff_file)
    path(genome_file)
    val(promotor_length)

    output:
    tuple val(meta), path("*_promotor.gff"), emit: promoter_bed
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python  /scratch/nadjafn/test_PROMOTOR_analysis/scripts/pormotor_extractor.py 
            --gff $gff_file \
            --fasta $genome_file \
            --output ${prefix}_promotor.gff \
            --length $promotor_length \
            --remove-overlaps

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gff_to_promotor: \$(python --version | sed -e "s/python v//g")
    END_VERSIONS
    """
}