
process MERGE_VCFS {
    tag "${meta.id}_${meta.synt_id}_${meta.direction}"
    publishDir "${params.outdir}/merged_vcfs", mode: 'copy'

    input:
    tuple val(meta), path(vcfs)

    output:
    tuple val(meta), path("${meta.id}_${meta.synt_id}_${meta.direction}_all_variants.vcf"), emit: merged_vcf
    path "versions.yml", emit: versions

    script:
    def prefix = "${meta.id}_${meta.synt_id}_${meta.direction}"
    """
    echo "##fileformat=VCFv4.2" > ${prefix}_all_variants.vcf
    grep --no-filename "##contig" ${vcfs} | sort | uniq >> ${prefix}_all_variants.vcf
    echo -e "#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\tSample" >> ${prefix}_all_variants.vcf
    grep --no-filename -v "#" ${vcfs} >> ${prefix}_all_variants.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | sed '1!d; s/.*version \\(.*\\)/\\1/')
    END_VERSIONS
    """
}