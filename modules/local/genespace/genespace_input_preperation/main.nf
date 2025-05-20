process GENESPACE_INPUT_PREPERATION {
    tag "$meta.id"
    label 'process_low'

    // conda "${moduleDir}/environment.yml"
    conda "/users/nadjafn/.conda/envs/orthofinder"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hdcf5f25_4' :
        'biocontainers/gffread:0.12.7--hdcf5f25_4' }"

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
        sed -i 's/\\(chr[0-9]\\+\\)_[0-9]\\+/\\1/g' ${prefix}/bed/\$(basename "\${gff%.*}").bed
    done

    for fasta in $fasta
    do
        cp \$fasta ${prefix}/peptide/\$(basename "\${fasta%.*}").fa
    done
        # remove the haplotype prefix

    """

}

