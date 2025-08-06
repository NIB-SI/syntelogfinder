include { GFF_TO_PROMOTOR } from '../../modules/local/gff_to_promotor'
include { PROMOTOR_EXTRACTION } from '../../modules/local/promotor_extraction'
include { MINIMAP2_ALIGN } from '../../modules/nf-core/minimap2/align/'
include { SVIM_ASM } from '../../modules/local/svim_asm/'
include { MERGE_VCFS} from '../../modules/local/merge_vcfs'
include { VCF_TO_SUMMARY } from '../../modules/local/vcf_to_summary'

workflow PROMOTOR_COMPARISON {
    take:
        gff_file
        fasta_ch
        promotor_length
        pangene_file
        synt_ids_ch

    main:
        promotor_gff = GFF_TO_PROMOTOR(gff_file, fasta_ch, promotor_length)

        // Combine the synt_ids with other inputs
        promotor_inputs = synt_ids_ch.flatten().combine(promotor_gff.gff)
                                     .combine(pangene_file)
                                     .combine(fasta_ch)

        promotors = PROMOTOR_EXTRACTION(promotor_inputs)

            // Generate all pairwise combinations
        pairwise_combinations = promotors.fasta
            .flatMap { id, synt_id, files ->
                def fasta_files = files.findAll { it.name.endsWith('.fasta') }
                def result = []
                for (int i = 0; i < fasta_files.size(); i++) {
                    for (int j = 0; j < fasta_files.size(); j++) {
                        if (i != j) {
                            result << [id, synt_id, fasta_files[i], fasta_files[j], 'query_to_reference']
                        }
                    }
                }
                return result
            }

        // Assuming your channel is named pairwise_combinations
        pairwise_combinations
            .map { id, synt_id, query, reference, direction ->
                def meta = [id: id, synt_id: synt_id, direction: direction, original_query: query.name, original_reference: reference.name]
                [meta, query, reference]
            }
            .set { minimap2_input }


        // Run minimap2 on pairwise combinations
        minimap2_pairwise_alignments = MINIMAP2_ALIGN(
            minimap2_input.map { it -> [it[0], it[1]] },  // [meta, reads]
            minimap2_input.map { it -> [it[0], it[2]] },
            true,  // bam_format
            "bai",   // bam_index_extension
            false,   // cigar_paf_format
            true   // cigar_bam
        )

        // Join the minimap2 output with the original pairwise_combinations
        svim_asm_input = minimap2_pairwise_alignments.bam.join(minimap2_pairwise_alignments.index, by:0).join(minimap2_input, by:0)


        svim_asm_output = SVIM_ASM(svim_asm_input)
        //svim_asm_output.vcf.groupTuple().view()

        collected_vcfs = svim_asm_output.vcf
            .map { meta, vcf ->
                def new_meta = [id: meta.id.id]
                [new_meta, meta.synt_id, vcf]
            }
            .groupTuple(by: [0,1])


        merged = MERGE_VCFS(
            collected_vcfs
        )


        summary_input = merged.merged_vcf
                .join(promotors.gff, by:[1])
                .map { synt_id, meta, vcf, meta2, gff ->
                        def new_meta = [id: meta]
                        [new_meta, synt_id, vcf, gff]
                        }

        VCF_TO_SUMMARY(
            summary_input,
            20
        )
        emit:
            summary = VCF_TO_SUMMARY.out.summary
}


