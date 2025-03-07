#!/usr/bin/env python3
import argparse
import pyranges as pr
import pandas as pd
import pyfaidx
import sys

def extend_features(gff_file, fasta_file, feature_type, extension, output_file):
    """
    Extend specified features in a GFF file and save the modified GFF.
    
    Args:
        gff_file (str): Path to input GFF file
        fasta_file (str): Path to reference FASTA file
        feature_type (str): Feature type to extend (e.g., 'CDS')
        extension (int): Number of base pairs to extend in both directions
        output_file (str): Path to output GFF file
    """
    # Read the GFF file
    print(f"Reading GFF file: {gff_file}")
    ann = pr.read_gff3(gff_file)
    
    # Select features of the specified type
    selector = (ann.Feature == feature_type)
    if not selector.any():
        print(f"No features of type '{feature_type}' found in the GFF file.")
        sys.exit(1)
    
    features = ann[selector]
    print(f"Found {len(features)} {feature_type} features")
    
    # Extend the features
    print(f"Extending {feature_type} features by {extension} bp")
    extended_features = features.extend({'5': extension, '3': extension}, group_by='Parent')
    
    # Remove original features and add extended ones
    ann_filtered = ann[~ann.Feature.isin([feature_type])]
    ann_updated_df = pd.concat([ann_filtered.df, extended_features.df])
    
    # Convert back to PyRanges and sort
    ann_updated = pr.PyRanges(ann_updated_df)
    ann_updated = ann_updated.sort(by='ID')
    
    # Clip features that extend beyond chromosome boundaries
    print(f"Clipping features to chromosome boundaries using: {fasta_file}")
    pyf = pyfaidx.Fasta(fasta_file)
    ann_updated = pr.genomicfeatures.genome_bounds(ann_updated, chromsizes=pyf, clip=True)
    
    # Write the updated GFF file
    print(f"Writing output to: {output_file}")
    ann_updated.to_gff3(output_file)
    print("Done!")

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Extend features in GFF file.')
    parser.add_argument('-g', '--gff', required=True, help='Input GFF file')
    parser.add_argument('-f', '--fasta', required=True, help='Reference FASTA file')
    parser.add_argument('-t', '--feature-type', default='CDS', help='Feature type to extend (default: CDS)')
    parser.add_argument('-e', '--extension', type=int, default=150, help='Extension length in bp (default: 150)')
    parser.add_argument('-o', '--output', required=True, help='Output GFF file')
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose output')
    
    args = parser.parse_args()
    
    # Run the feature extension
    try:
        extend_features(args.gff, args.fasta, args.feature_type, args.extension, args.output)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()




# ./extend_gff_features.py \
#   --gff /scratch/nadjafn/potato-allelic-orthogroups/results/0_LongestIsofrom/De_v1.functional_annotations_nodigits.longest_isoforms.gff \
#   --fasta "/scratch/nadjafn/reference/Desiree_v1/De_v1.fa" \
#   --feature-type CDS \
#   --extension 150 \
#   --output /scratch/nadjafn/reference/Desiree_v1/De_v1.functional_annotations_nodigits.longest_isoforms_cds150.gff