// modules/gff_to_promoter.nf

process GFF_TO_PROMOTOR {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"


    input:
    tuple val(meta), path(gff_file)
    path(genome_file)
    val(promotor_length)

    output:
    tuple val(meta), path("*_promotor.gff"), emit: gff
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python ${baseDir}/scripts/promotor_extractor.py \
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
