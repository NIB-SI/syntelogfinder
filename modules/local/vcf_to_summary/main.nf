process VCF_TO_SUMMARY {
    tag "$meta.id"
    label 'process_medium'

    conda "/users/nadjafn/.conda/envs/syri"
    //conda "bioconda::pysam=0.19.1 conda-forge::matplotlib=3.7.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:6684529f3d3b42b8b428edeb1d3e5c5c3bce1d7d-0' :
        'biocontainers/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:6684529f3d3b42b8b428edeb1d3e5c5c3bce1d7d-0' }"

    input:
    tuple val(meta), val(synt_id), path(vcf), path(gff)
    val min_sv_length

    output:
    tuple val(meta), val(synt_id), path("${prefix}*synt_id_summary.csv"), emit: summary
    tuple val(meta), val(synt_id), path("${prefix}*.svg"), emit: plots, optional: true
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_${synt_id}"
    """
    # merge gff files
    cat ${gff} | grep -v "#" | sort -k1,1 -k4,4n > ${prefix}_promotor.gff
    python /scratch/nadjafn/test_PROMOTOR_analysis/scripts/vcf_to_summary.py \\
        --vcf $vcf \\
        --gff  ${prefix}_promotor.gff \\
        --min-sv-length $min_sv_length \\
        --plot \\
        --out ${prefix} \\
        --synt_id ${synt_id} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pysam: \$(python -c "import pysam; print(pysam.__version__)")
        matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """
}