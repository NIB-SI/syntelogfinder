process SPLIT_HAPLOTYPES {
    tag "$meta.id"
    label 'process_single'

    // conda "${moduleDir}/environment.yml"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/agat:1.4.0--pl5321hdfd78af_0' :
    //     'biocontainers/agat:1.4.0--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("hap*.gff"), emit: output_gtf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ploidy = task.ext.ploidy ?: 4
    """
    for ((i=1; i<=$ploidy; i++)); do
        grep -e "chr_\\?[0-9]\\+_\${i}" "$gff" > "hap\${i}.gff"
        # remove the _haplotype suffix on chr name
        # sed 's/\\(chr[0-9]\\+\\)_[0-9]\\+/\\1/g' "hap\${i}_pre.gff" > "hap\${i}.gff"
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