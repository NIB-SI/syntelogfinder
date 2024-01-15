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
include { run_diamond_all } from './Modules/Diamond'
include { run_MCScanX } from './Modules/MCscanX'


params.reference_dir = '/scratch/nadjafn/potato-allelic-orthogroups/potato_culitvar_references'
params.mcScanX_dir = 'Input_MCscanX/cultivars'



// genome_with_gff = Channel.
//     fromFilePairs("${params.reference_dir}/*{.asm.fa,.hc_gene_models.repr.gff3}")
//     | map { id, reads ->
//     [id, reads[0], reads[1]]
//     }   | view


prot_dir = Channel.fromPath(params.prot_dir)
mcScanX_dir = Channel.fromPath(params.mcScanX_dir)
culitvar = Channel.of('AT_CR_OT')


workflow {
    //get_CDS(genome_with_gff)
    //run_OrthoFinder(prot_dir)
    (gff_dir, prot_dir) = download_rename(culitvar)
    diamond_blast = run_diamond_all(prot_dir)
    merged_bed = merge_gff_to_bed(gff_dir)
    run_MCScanX(merged_bed, diamond_blast)
    run_OrthoFinder(prot_dir)
}
