/* 
 * pipeline input parameters 
 */


log.info """\
         ASE benchmarking P I P E L I N E    
         ===================================
         
         """
         .stripIndent()

// import processes from modules
include { download_rename; merge_gff_to_bed } from './Modules/PrepareProteoms'
include { run_OrthoFinder } from './Modules/OrthoFinder'
include { run_diamond_all; filter_diamond } from './Modules/Diamond'
include { run_MCScanX } from './Modules/MCscanX'
include { run_cd_hit } from './Modules/CDHit'
include { add_parent_to_gff; run_tranD; gff3_to_gtf; count_exons } from './Modules/TranD'
include { run_blastn } from './Modules/Blast'
include { look_at_exact_dups } from './Modules/Analysis'


params.reference_dir = '/scratch/nadjafn/potato-allelic-orthogroups/potato_culitvar_references'
params.mcScanX_dir = 'Input_MCscanX/cultivars'
params.hc_gene_models_gff = '/scratch/nadjafn/potato-allelic-orthogroups/gff_reference/ATL_v3.hc_gene_models.gff3'



// genome_with_gff = Channel.
//     fromFilePairs("${params.reference_dir}/*{.asm.fa,.hc_gene_models.repr.gff3}")
//     | map { id, reads ->
//     [id, reads[0], reads[1]]
//     }   | view


prot_dir = Channel.fromPath(params.prot_dir)
mcScanX_dir = Channel.fromPath(params.mcScanX_dir)
culitvar = Channel.of('AT_CR_OT')
hc_gene_models_gff  = Channel.fromPath(params.hc_gene_models_gff)

workflow {
    //get_CDS(genome_with_gff)
    (gff_dir, prot_dir, cDNA_dir) = download_rename(culitvar)
    diamond_blast = run_diamond_all(prot_dir)

    merged_bed = merge_gff_to_bed(gff_dir)
    run_MCScanX(merged_bed, diamond_blast)
    run_OrthoFinder(prot_dir, 'A')
    diamond_out = filter_diamond(prot_dir, 'A')
    blast_out = run_blastn(cDNA_dir, 'A')

    // analyse number of dublicated genes
    look_at_exact_dups(blast_out.concat(diamond_out))

    gff = add_parent_to_gff(hc_gene_models_gff)
    gtf_file = gff3_to_gtf(gff)
    gtf_file = count_exons(gtf_file)
    //run_tranD(gtf_file)

}
