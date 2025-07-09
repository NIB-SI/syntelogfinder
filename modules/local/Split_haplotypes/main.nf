process SPLIT_HAPLOTYPES {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("hap*.gff"), emit: output_gtf
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Extract unique haplotype suffixes from the GFF file
    # This works for patterns like BrgdChr01A, BrgdChr01B, etc.
    haplotypes=\$(grep -o 'BrgdChr[0-9]\\+[A-Z]' "$gff" | \\
                 sed 's/.*\\([A-Z]\\)\$/\\1/' | \\
                 sort -u)

    # If no BrgdChr pattern found, try the original numeric pattern
    if [[ -z "\$haplotypes" ]]; then
        haplotypes=\$(seq 1 $params.ploidy)
    fi
    

    # Split into haplotype-specific files
    for hap in \$haplotypes; do
        if [[ "\$hap" =~ ^[0-9]+\$ ]]; then
            # Numeric haplotype (original pattern)
            grep -e "chr_\\?[0-9]\\+_\${hap}" "$gff" > "hap\${hap}.gff"
        else
            # Letter-based haplotype (new pattern)
            grep "BrgdChr[0-9]\\+\${hap}" "$gff" > "hap\${hap}.gff"
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.agat.gtf
    touch ${gff}.agat.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    agat: \$(agat_sp_keep_longest_isoform.pl --help | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p')
    END_VERSIONS
    """
}
