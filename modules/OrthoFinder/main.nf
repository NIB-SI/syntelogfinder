process run_OrthoFinder {
    publishDir(
        path: "${params.publish_dir}/OrthoFinder",
    )
    cache 'lenient'
    conda "/users/nadjafn/.conda/envs/orthofinder"
    input:
        tuple val(meta), path(pep)

    output:
        path("Orthofinder_output")
    script:
    """
        # split pep into haplotypes 
        awk '
    /^>/ {
        # Check the header line for H1, H2, H3, or H4.
        if (\$0 ~ "H1") {
        out="H1_sequences.fasta"
        } else if (\$0 ~ "H2") {
        out="H2_sequences.fasta"
        } else if (\$0 ~ "H3") {
        out="H3_sequences.fasta"
        } else if (\$0 ~ "H4") {
        out="H4_sequences.fasta"
        } 
    print \$0 > out
    next
    }
    {
        # For sequence lines (non-header), just continue printing to the current file.
        print \$0 > out
    }
    ' input.fasta
    mkdir -p filtered_prot_dir
    cp *sequences.fasta filtered_prot_dir 
    orthofinder -t ${task.cpu} -f filtered_prot_dir -o Orthofinder_output
    """
}