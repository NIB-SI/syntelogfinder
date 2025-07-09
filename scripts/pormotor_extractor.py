#!/usr/bin/env python3
"""
Promoter Region Extractor

This script extracts promoter regions from genomic data based on GFF annotations.
It can optionally remove overlaps with coding sequences.

Usage:
    python promoter_extractor.py --gff <gff_file> --fasta <fasta_file> \
                                 --output <output_file> --length <promoter_length> \
                                 [--remove-overlaps]
"""

import argparse
import os
import sys
import pyranges as pr
import pandas as pd
import pyfaidx
from pathlib import Path


def extract_promoters(gff_file, fasta_file, output_file, promoter_length=3000, remove_overlaps=False):
    """
    Extract promoter regions from GFF file based on CDS features.
    
    Args:
        gff_file (str): Path to the GFF file
        fasta_file (str): Path to the FASTA file
        output_file (str): Path to the output GFF file
        promoter_length (int): Length of the promoter region to extract
        remove_overlaps (bool): Whether to remove overlaps with coding sequences
    """
    print(f"Reading GFF file: {gff_file}")
    ann = pr.read_gff3(gff_file)
    
    # Select CDS features
    print(f"Extracting CDS features")
    selector = (ann.Feature == 'CDS')
    cds = ann[selector]
    
    # Extend CDS features to include promoter regions
    print(f"Extending features by {promoter_length}bp")
    g = cds.extend({'5': promoter_length}, group_by='Parent')
    g.Feature = 'exon'
    
    # Extract promoter regions
    print(f"Extracting promoter sequences")
    prom = g.spliced_subsequence(0, promoter_length, "Parent")
    
    # Correct bounds to chromosome size
    print(f"Correcting bounds to chromosome size")
    pyf = pyfaidx.Fasta(fasta_file)
    prom = pr.genomicfeatures.genome_bounds(prom, chromsizes=pyf, clip=True)
    
    if remove_overlaps:
        print(f"Removing overlaps with coding sequences")
        # Check the overlap of promoters and CDS
        prom_in_cds = prom.intersect(cds, strandedness=False)
        
        # If there are overlaps, adjust the promoters
        if not prom_in_cds.df.empty:
            prom_in_cds_merged = prom_in_cds.merge(by='Parent', slack=10000)
            prom_in_cds_merged_df = prom_in_cds_merged.df
            
            # For positive strand promoters, set the start to be the end of the overlap
            mask_positive = prom_in_cds_merged_df['Strand'] == '+'
            prom_in_cds_merged_df.loc[mask_positive, 'Start'] = prom_in_cds_merged_df.loc[mask_positive, 'End']
            
            # For negative strand promoters, set the end to be the start of the overlap
            mask_negative = prom_in_cds_merged_df['Strand'] == '-'
            prom_in_cds_merged_df.loc[mask_negative, 'End'] = prom_in_cds_merged_df.loc[mask_negative, 'Start']
            
            # Get the positive strand features
            prom_in_cds_merged_pos = prom_in_cds_merged_df[prom_in_cds_merged_df.Strand == '+']
            prom_overlap_start_pos = prom_in_cds_merged_pos[["Start", "Parent"]]
            
            # Get the negative strand features
            prom_in_cds_merged_neg = prom_in_cds_merged_df[prom_in_cds_merged_df.Strand == '-']
            prom_overlap_stop_neg = prom_in_cds_merged_neg[["End", "Parent"]]
            
            # Update the promoter dataframe
            prom_df = prom.df
            prom_df.set_index('Parent', inplace=True)
            prom_overlap_start_pos.set_index('Parent', inplace=True)
            prom_overlap_stop_neg.set_index('Parent', inplace=True)
            
            prom_df.update(prom_overlap_start_pos, overwrite=True)
            prom_df.update(prom_overlap_stop_neg, overwrite=True)
            
            # Reset index and prepare for output
            prom_df["Parent"] = prom_df.index
            prom_df.reset_index(drop=True, inplace=True)
        else:
            prom_df = prom.df
    else:
        prom_df = prom.df
    
    # Sort and prepare for output
    prom_df.sort_values(by=['Chromosome', 'Start'], inplace=True)
    
    # Add Name and ID fields for the promoters
    if 'Parent' in prom_df.columns:
        prom_df["Name"] = prom_df["Parent"]
        prom_df["ID"] = "promoter_" + prom_df["Parent"]
        # Drop Parent column to avoid problems with relationships
        prom_df.drop(columns=['Parent'], inplace=True)
    else:
        print("Warning: No 'Parent' column found in the dataframe.")
    
    # Create PyRanges object and set feature type
    prom_updated = pr.PyRanges(prom_df)
    prom_updated.Feature = 'exon'
    
    # Write to output file
    print(f"Writing output to: {output_file}")
    prom_updated.to_gff3(output_file)
    
    # Print summary statistics
    print(f"\nSummary:")
    print(f"  Total promoters extracted: {len(prom_df)}")
    if remove_overlaps:
        print(f"  Promoters adjusted for CDS overlaps: {len(prom_in_cds_merged_df) if 'prom_in_cds_merged_df' in locals() else 0}")
    
    return prom_updated


def main():
    parser = argparse.ArgumentParser(description='Extract promoter regions from genomic data.')
    parser.add_argument('--gff', required=True, help='Input GFF file')
    parser.add_argument('--fasta', required=True, help='Input FASTA file')
    parser.add_argument('--output', required=True, help='Output GFF file')
    parser.add_argument('--length', type=int, default=3000, help='Promoter length (default: 3000bp)')
    parser.add_argument('--remove-overlaps', action='store_true', help='Remove overlaps with coding sequences')
    
    args = parser.parse_args()
    
    # Check if input files exist
    if not Path(args.gff).is_file():
        sys.exit(f"Error: GFF file '{args.gff}' not found.")
    if not Path(args.fasta).is_file():
        sys.exit(f"Error: FASTA file '{args.fasta}' not found.")
    
    # Check if output directory exists
    output_dir = Path(args.output).parent
    if not output_dir.exists():
        sys.exit(f"Error: Output directory '{output_dir}' does not exist.")
    
    # Extract promoters
    try:
        extract_promoters(
            args.gff,
            args.fasta,
            args.output,
            args.length,
            args.remove_overlaps
        )
        print("Promoter extraction completed successfully.")
    except Exception as e:
        sys.exit(f"Error: {str(e)}")


if __name__ == "__main__":
    main()



