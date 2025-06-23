process PROMOTOR_EXTRACTION {
    tag "$meta.id"
    label 'process_low'

    conda "/users/nadjafn/.conda/envs/gffread"


    input:
    tuple val(synt_id), val(meta), path(promoter_gff_file), val(meta2), val(haplotypes), path(synt_file), path(genome_file)

    output:
    tuple val(meta), val(synt_id), path("individual_seqs/*fasta"), emit: fasta
    tuple val(meta.id), val(synt_id), path("individual_seqs/*gff"), emit: gff
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p individual_seqs

    # Get the gene IDs from the syntelog file
    gene_ids=\$(grep -w ${synt_id} ${synt_file} | awk '{print \$1}')

    # Process each gene ID
    for gene_id in \$gene_ids; do
        grep "promoter_\${gene_id}" ${promoter_gff_file} > individual_seqs/\${gene_id}_promotor.gff
        gffread -g ${genome_file} -w individual_seqs/\${gene_id}_promotor.fasta individual_seqs/\${gene_id}_promotor.gff
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1 | sed 's/^.*gffread v//; s/ .*\$//')
    END_VERSIONS
    """
}
