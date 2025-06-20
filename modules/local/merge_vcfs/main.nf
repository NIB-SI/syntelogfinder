
process MERGE_VCFS {
    tag "${meta.id}_${synt_id}"
    publishDir "${params.outdir}/merged_vcfs", mode: 'copy'

    input:
    tuple val(meta), val(synt_id), path(vcfs)

    output:
    tuple val(meta.id), val(synt_id), path("${meta.id}_${synt_id}_all_variants.vcf"), emit: merged_vcf
    path "versions.yml", emit: versions

    script:
    def prefix = "${meta.id}_${synt_id}"
    def vcf_file_list = vcfs.collect { it.toString() }.join(' ')
    """
    # Check if all input files exist
    for vcf in $vcf_file_list; do
        if [ ! -f "\$vcf" ]; then
            echo "Error: Input file \$vcf does not exist"
            exit 1
        fi
    done

    echo "##fileformat=VCFv4.2" > ${prefix}_all_variants.vcf
    grep --no-filename "##contig" $vcf_file_list | sort | uniq >> ${prefix}_all_variants.vcf
    echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample" >> ${prefix}_all_variants.vcf
    grep --no-messages --no-filename -v "#" $vcf_file_list >> ${prefix}_all_variants.vcf || true
    
    cat <<-END_VERSIONS > versions.yml
    "PROMOTOR_COMPARISON:MERGE_VCFS":
        bash: \$(bash --version | sed '1!d; s/.*version \\(.*\\)/\\1/')
    END_VERSIONS
    """
}