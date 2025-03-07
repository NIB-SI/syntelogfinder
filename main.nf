#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Log pipeline information
log.info """
         Syntleog P I P E L I N E    
         ===================================
         """
         .stripIndent()

// Import processes from modules
include { AGAT_spKeepLongestIsoform } from './modules/local/AGAT/spKeepLongestIsoform'
include { GFFREAD as GFFREAD_PROT; GFFREAD as GFFREAD_BED } from './modules/local/gffread'
include { SPLIT_HAPLOTYPES } from './modules/local/Split_haplotypes'
include { GENESPACE_INPUT_PREPERATION } from './modules/local/genespace/genespace_input_preperation'
include { GENESPACE_RUN } from './modules/local/genespace/genespace_run'
include { GENESPACE_PARSE } from './modules/local/genespace/genespace_parse'
include { EXTEND_GFF_FEATURES } from './modules/local/extend_gff_features'
include { CDS_BLAST } from './subworkflows/cds_blast'

// Define pipeline parameters

    // Input files
params.reference_fasta = '/scratch/nadjafn/reference/Desiree_v1/De_v1.fa'
params.reference_gff   = '/scratch/nadjafn/reference/Desiree_v1/De_v1.functional_annotations_nodigits.gff'
    
    // Output directory
params.outdir = '/DKED/scratch/nadjafn/potato-allelic-orthogroups/output'
    
    // Tool paths
params.mcscanx_path = '/DKED/scratch/nadjafn/MCScanX'
    
    // Config file
params.config = '/DKED/scratch/nadjafn/potato-allelic-orthogroups/conf/nextflow.config'


// Create input channels
gff_ch = Channel.fromPath(params.reference_gff)
    .map { gff -> tuple([id: gff.baseName], gff) }

fasta_ch = Channel.fromPath(params.reference_fasta)

// Main workflow
workflow {
    // Keep longest isoform
    agat_output = AGAT_spKeepLongestIsoform(gff_ch)

    // Split GFF into haplotypes
    split_gff = SPLIT_HAPLOTYPES(agat_output.output_gtf)

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

    genespace_input = GENESPACE_INPUT_PREPERATION(gffread_output)

    // Run GENESPACE
    genespace_run = GENESPACE_RUN(genespace_input.dir, params.mcscanx_path)

    // Parse GENESPACE output
    genespace_parse = GENESPACE_PARSE(genespace_run.pangenes.join(agat_output.output_gtf))

    // Extend GFF features
    extended_gff = EXTEND_GFF_FEATURES(agat_output.output_gtf, fasta_ch)

    // Run CDS BLAST subworkflow
    CDS_BLAST(extended_gff.gff, haplotype_ch.map{it.fasta})
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