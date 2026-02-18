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
    # Extract unique haplotype suffixes from the GFF file.
    haplotypes=\$(grep -o '[Cc]hr_\\?[0-9]\\+_\\?[A-Z1-9]' "$gff" | \
    sed 's/.*\\([A-Z1-9]\\)\$/\\1/' | \
    sort -u)

    # If no Chr pattern found, try the original numeric pattern
    if [[ -z "\$haplotypes" ]]; then
        haplotypes=\$(seq 1 $params.ploidy)
    fi

    # Split into haplotype-specific files
    for hap in \$haplotypes; do
        echo \$hap
        if [[ "\$hap" =~ ^[1-9]+\$ ]]; then
            # Numeric haplotype (original pattern)
            grep "chr_\\?[0-9]\\+_\${hap}\\b" "$gff" > "hap\${hap}.gff"
        else
            # Alphanumeric haplotype (new pattern)
            grep "[Cc]hr[0-9]\\+_\\?\${hap}\\b" "$gff" > "hap\${hap}.gff"
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        grep: \$(grep --version | head -n1 | sed 's/.*[[:space:]]\\([0-9][0-9.]*\\).*/\\1/')
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
        grep: \$(grep --version | head -n1 | sed 's/.*[[:space:]]\\([0-9][0-9.]*\\).*/\\1/')
    END_VERSIONS
    """
}
