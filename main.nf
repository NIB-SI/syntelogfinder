#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Log pipeline information
log.info """
         Syntelog P I P E L I N E
         ===================================
         """
         .stripIndent()

// Import processes from modules
include { AGAT_spKeepLongestIsoform } from './modules/local/AGAT/spKeepLongestIsoform'
include { GFFREAD as GFFREAD_PROT; GFFREAD as GFFREAD_BED } from './modules/nf-core/gffread'
include { SPLIT_HAPLOTYPES } from './modules/local/Split_haplotypes'
include { GENESPACE_ANALYSIS } from './subworkflows/genespace_analysis'
include { PROMOTOR_COMPARISON } from './subworkflows/promotor_comparison'
include { EXTEND_GFF_FEATURES } from './modules/local/extend_gff_features'
include { CDS_BLAST } from './subworkflows/cds_blast'
include { SYNTELOG_SIMILARITY } from './modules/local/Syntelog_similarity'

// Define pipeline parameters


    // Tool paths
params.mcscanx_path = '/DKED/scratch/nadjafn/MCScanX'

    // Config file
params.config = "${baseDir}/conf/nextflow.config"


// Create input channels
gff_ch = Channel.fromPath(params.reference_gff)
    .map { gff -> tuple([id: gff.baseName], gff) }

fasta_ch = Channel.fromPath(params.reference_fasta)

promotor_length = Channel.from([3000])

// Main workflow
workflow {
    // Keep longest isoform
    agat_output = AGAT_spKeepLongestIsoform(gff_ch)

    // Split GFF into haplotypes
    split_gff = SPLIT_HAPLOTYPES(agat_output.output_gff)

    // Prepare haplotype channels
    haplotype_ch = split_gff.output_gtf
        .transpose()
        .combine(fasta_ch)
        .map { meta, gff, fasta ->
            [
                gff: tuple([id: gff.baseName, cultivar: meta], gff),
                fasta: fasta
            ]
        }

    // Extract protein sequences and convert GFF to BED
    gffread_output_prot = GFFREAD_PROT(haplotype_ch.map{it.gff}, haplotype_ch.map{it.fasta})
    gffread_output_bed = GFFREAD_BED(haplotype_ch.map{it.gff}, haplotype_ch.map{it.fasta})


    // Prepare GENESPACE input
    gffread_output = gffread_output_prot.gffread_fasta
        .map { meta, gff -> tuple(meta.cultivar, gff.baseName, gff) }
        .join(
            gffread_output_bed.gffread_gff
                .map { meta, gff -> tuple(meta.cultivar, gff.baseName, gff) },
            by: [0, 1]
        )
        .groupTuple(by: 0)

    // GENESPACE Analysis
    genespace_ch = GENESPACE_ANALYSIS(gffread_output, agat_output.output_gff)

    // Extend GFF features
    extended_gff = EXTEND_GFF_FEATURES(gff_ch, fasta_ch)

    if (params.run_blast) {
        // Check if the CDS_BLAST subworkflow is enabled
        log.info "Running CDS_BLAST subworkflow"
        // Run CDS BLAST subworkflow
        blast_ch = CDS_BLAST(extended_gff.gff, haplotype_ch.map{it.fasta})


        SYNTELOG_SIMILARITY(
        genespace_ch,
        blast_ch.results)

    } else {
        log.info "Skipping CDS_BLAST subworkflow"
    }

    if (params.run_promotor_comparison) {
        // Check if the Syntelog similarity module is enabled
        log.info "Running Syntelog similarity module"
        // Run Syntelog similarity module
        genespace_ch.view()
        // Define the synteny category you want to filter for
        params.synteny_category = "1hap1_1hap2_1hap3_1hap4_s"
        // Process the file and create a new channel with filtered Synt_ids
        synt_id_ch = genespace_ch
            .map { meta, haplotypes, file -> file }  // Extract just the file path
            .splitCsv(sep: '\t', header: true)
            .filter { row -> row.synteny_category == params.synteny_category }
            .map { row -> row.Synt_id }
            .unique()
            .take(10)  // Take only the first 10 unique Synt_ids
            .collect()


        // View the results (for debugging)
        synt_id_ch.view()
        // Run Promotor comparision subworkflow
        synt_id_ch = Channel.from(['Synt_id_15856','Synt_id_17129'])

        promotor_summary = PROMOTOR_COMPARISON(agat_output.output_gff, fasta_ch, promotor_length, genespace_ch, synt_id_ch)
        promotor_summary.view()
    } else {
        log.info "Skipping Promotor comparision subworkflow"
    }

}

// Function to check if required parameters are set
def checkRequiredParams() {
    def required = ['reference_fasta', 'reference_gff', 'outdir', 'mcscanx_path']
    for (param in required) {
        if (params[param] == null) {
            exit 1, "Required parameter '${param}' is missing"
        }
    }
}

// Call the function to check required parameters
checkRequiredParams()
