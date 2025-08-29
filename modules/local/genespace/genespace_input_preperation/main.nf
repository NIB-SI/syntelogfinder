process GENESPACE_INPUT_PREPERATION {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), val(haplotypes), path(fasta), path(gff)

    output:
    tuple val(meta), val(haplotypes), path(fasta), path(gff), path("${meta.id}_genespace"), emit: dir             , optional: false

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args             ?: ''
    def prefix      = task.ext.prefix           ?: "${meta.id}_genespace"
    """
    mkdir -p ${prefix}/bed
    mkdir -p ${prefix}/peptide
    for gff in $gff
    do
        grep "CDS" \$gff | cut -f1-4 > ${prefix}/bed/\$(basename "\${gff%.*}").bed
        # remove haplotype suffix to match the same chromosome on different haplotypes for GENESPACE sameChr option
        sed -i 's/^\\([A-Za-z]*[Cc]hr_\\?[0-9]\\{1,2\\}\\)[^\\t]*/\\1/'  ${prefix}/bed/\$(basename "\${gff%.*}").bed
    done

    for fasta in $fasta
    do
        cp \$fasta ${prefix}/peptide/\$(basename "\${fasta%.*}").fa
        # remove the haplotype prefix to match gff file

    done


    """

}

