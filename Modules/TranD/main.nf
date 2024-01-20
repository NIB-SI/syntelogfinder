process add_parent_to_gff {

    cache 'lenient'
    input:
        path(gff_file)
    output:
        path("with_parent.gff3")
    script:
    """
    # Process the file
    awk '{
        if(\$3 == "mRNA") {
            # Find the Name field and extract the ID without the last .number
            match(\$0, /Name=[^;]+/);
            name_field = substr(\$0, RSTART, RLENGTH);
            split(name_field, name_arr, "=");
            split(name_arr[2], id_arr, "\\.");
            id = id_arr[1];
            for(i = 2; i < length(id_arr); i++) {
                id = id"."id_arr[i];
            }
            # Add the Parent field
            print \$0";Parent="id;
        } else {
            print \$0;
        }
    }' ${gff_file} > with_parent.gff3
    """
}


process gff3_to_gtf {

    cache 'lenient'
    input:
        path(gff_file)
    conda "/users/nadjafn/.conda/envs/gffread"

    output:
        path("my.gtf")
    script:
    """
    gffread ${gff_file} --keep-genes -F -T -o my.gtf
    """
}

process count_exons {
    cache 'lenient'
    input:
        path(gtf_file)
    conda "/users/nadjafn/.conda/envs/pygtftk"

    output:
        tuple path("with_genes.gtf"),
              path("*_count.gtf")
    script:
    """
    gtftk convert_ensembl -i $gtf_file > with_genes.gtf
    gtftk nb_exons -i with_genes.gtf | gtftk nb_transcripts -o nb_transcripts.out.gtf
    gtftk tabulate -k gene_id,seqid,start,end,nb_tx  -i nb_transcripts.out.gtf > gene_count.gtf
    gtftk tabulate -k transcript_id,seqid,start,end,nb_exons  -i nb_transcripts.out.gtf > transcript_count.gtf
    """
}
// cat ${gff_dir}/ATL_v3.working_models.gff3 | grep 'exon' | awk -F'\\t' '{
//     split(\$9, a, ";");
//     parent_id = "";
//     for (i in a) {
//         if (match(a[i], /Parent=/)) {
//             parent_id = substr(a[i], RSTART + 7);
//         }
//     }
//     n = split(parent_id, b, ".");
//     gene_id = b[1];
//     for (i = 2; i <= 3 && i <= n; i++) {
//         gene_id = gene_id "." b[i];
//     }
//     \$9 = "gene_id \\"" gene_id "\\"; transcript_id \\"" parent_id "\\"";
//     print
// }' OFS='\\t' | sed 's/\\tParent=[^\\t]*//' > converted.gff3


process run_tranD {
    publishDir(
        path: "${params.publish_dir}/tranD",
    )
    cache 'lenient'
    conda "/users/nadjafn/.conda/envs/tranD_env"
    input:
        path(gtf_file)

    output:
        path("tranD_results")
    script:
    """
    trand $gtf_file --outdir tranD_results --cpus 16
    """
}


// I think a file with the gene name start end number of exons , number of  transcripts will be very useful.