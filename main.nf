/* 
 * pipeline input parameters 
 */


log.info """\
         Syntleog P I P E L I N E    
         ===================================
         
         """
         .stripIndent()

// import processes from modules
include { download_rename; merge_gff_to_bed } from './modules/PrepareProteoms'
include { run_OrthoFinder } from './modules/OrthoFinder'
include { DIAMOND } from './modules/Diamond'
include { run_MCScanX } from './modules/MCscanX'
include { run_cd_hit } from './modules/CDHit'
include { run_blastn } from './modules/Blast'
include { look_at_exact_dups } from './modules/Analysis'
include { AGAT_spKeepLongestIsoform} from './modules/AGAT/spKeepLongestIsoform'
include { GFFREAD as GFFREAD_PROT; GFFREAD as GFFREAD_BED}  from './modules/gffread'
include { SPLIT_HAPLOTYPES } from './modules/Split_haplotypes'
include { GENESPACE_INPUT_PREPERATION } from './modules/genespace/genespace_input_preperation'
include { GENESPACE_RUN } from './modules/genespace/genespace_run'
include { GENESPACE_PARSE } from './modules/genespace/genespace_parse'
include { CDS_BLAST } from './subworkflows/cds_blast'
// problem:: task.ext. only applicable to GFFREAD
// maybe use bed tools to convert gff to bed
params.reference_fasta = '/scratch/nadjafn/reference/Desiree_v1/De_v1.fa'
// params.reference_fasta = "/scratch/nadjafn/reference/Atlantic/ATL_v3.asm.fa"
params.reference_gff = '/scratch/nadjafn/reference/Desiree_v1/De_v1.functional_annotations_nodigits.gff'
// params.reference_gff = "/scratch/nadjafn/reference/Atlantic/ATL_v3.hc_gene_models.repr.gff3"
// params.reference_gff = "/scratch/nadjafn/reference/Atlantic/updated_v62Atl_liftoff_a0.9s0.9_ALL.gff"
params.outdir = '/DKED/scratch/nadjafn/potato-allelic-orthogroups/output'
params.mcscanx_path = '/DKED/scratch/nadjafn/MCScanX'
params.config = '/DKED/scratch/nadjafn/potato-allelic-orthogroups/conf/nextflow.config'




// genome_with_gff = Channel.
//     fromFilePairs("${params.reference_dir}/*{.asm.fa,.hc_gene_models.repr.gff3}")
//     | map { id, reads ->
//     [id, reads[0], reads[1]]
//     }   | view


gff  = Channel.fromPath(params.reference_gff)


gff.map { gff ->
    tuple(meta = [id: gff.baseName], gff)
}.set { gff_ch }

fasta = Channel.fromPath(params.reference_fasta)

workflow {
    //get_CDS(genome_with_gff)

    // run AGAT_spKeepLongestIsoform
    agat_output = AGAT_spKeepLongestIsoform(gff_ch)
    agat_output.output_gtf.view()
    // split gff into haplotypes
    split_gff = SPLIT_HAPLOTYPES(agat_output.output_gtf)

    // put haplotypes in seperate channels
    haplotype_ch = split_gff.output_gtf.transpose().combine(fasta)
    haplotype_ch.multiMap { id, gff, fasta ->
        gff: tuple(meta = [id: gff.baseName, cultivar: id], gff)
        fasta: [fasta]
    }.set { haplotype_ch }

    // run gffread to extract protein sequences
    gffread_output_prot = GFFREAD_PROT(haplotype_ch.gff, haplotype_ch.fasta)

    // run gffread to convert gff to bed 
    gffread_output_bed = GFFREAD_BED(haplotype_ch.gff, haplotype_ch.fasta)

    // combine gffread outputs for genespace
    gffread_output_prot.gffread_fasta.map { meta, gff ->
        tuple(meta = meta.cultivar, gff.baseName, gff)
    }.set { gffread_output_prot }


    gffread_output_bed.gffread_gff.map { meta, gff ->
        tuple(meta = meta.cultivar, gff.baseName, gff)
    }.set { gffread_output_bed }

    // Combine the two channels
    gffread_output = gffread_output_prot.join(gffread_output_bed, by: [0, 1]).groupTuple(by: 0)

    // prepare genespace run
    genespace_input = GENESPACE_INPUT_PREPERATION(gffread_output)

    // run genespace
    genespace_run = GENESPACE_RUN(genespace_input.dir, params.mcscanx_path)

    // parse genespace output
    genespace_parse = GENESPACE_PARSE(genespace_run.pangenes.join(agat_output.output_gtf), 300)

    // Run the subworkflow
    CDS_BLAST(
        genespace_parse.plusgff,
        haplotype_ch.fasta
    )

}
