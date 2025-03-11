#!/usr/bin/env python3
"""
Syntenic Gene Analysis Tool

This script analyzes BLAST results and syntenic genes to identify allelic relationships
between haplotypes. It categorizes syntenic groups based on sequence identity patterns
and visualizes the distribution of SNPs and allelic categories.
"""

import pandas as pd
import argparse
import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple
from pathlib import Path


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Analyze BLAST results and syntenic genes')
    parser.add_argument('-b', '--blast', help='BLAST all-vs-all output file (fmt6)', required=True)
    parser.add_argument('-s', '--syntenic_genes', help='TSV file with syntenic genes', required=True)
    parser.add_argument('-o', '--output_prefix', help='Prefix for output files', required=True)
    
    return parser.parse_args()


def parse_blast_fmt6(blast_file: str) -> pd.DataFrame:
    """
    Parse BLAST output file in tabular format 6 into a pandas DataFrame.
    
    Args:
        blast_file: Path to BLAST output file
        
    Returns:
        DataFrame with parsed BLAST results
    """
    columns = [
        'query', 'subject', 'identity', 'length', 'mismatch', 'gap',
        'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bitscore'
    ]
    
    df = pd.read_csv(blast_file, sep='\t', header=None, names=columns)
    
    # Remove self-matches (same query and subject)
    df = df[df['query'] != df['subject']]
    # remove duplicated query subject pairs
    df = df.drop_duplicates(subset=['query', 'subject'])
    
    return df


def load_syntelogs(syntelogs_file: str) -> pd.DataFrame:
    """
    Load and prepare syntelogs data for analysis.
    
    Args:
        syntelogs_file: Path to syntelogs TSV file
        
    Returns:
        Processed syntelogs DataFrame with equal-length entries only
    """
    # Load syntelogs data
    syntelogs = pd.read_csv(
        syntelogs_file,
        sep='\t',
        usecols=['transcript_id', 'Synt_id', 'synteny_category', 'CDS_haplotype_with_longest_annotation', 'haplotype', 'CDS_length_category']
    )
    return syntelogs

def filter_syntelogs(syntelogs: pd.DataFrame) -> pd.DataFrame:
    """
    Filter syntelogs for equal lengths and extract haplotype information.
    """ 
    # Filter for equal lengths and extract haplotype information
    syntelogs = syntelogs[syntelogs['CDS_haplotype_with_longest_annotation'] == 'equal_lengths']
    # drop duplicated transcript_id
    syntelogs = syntelogs.drop_duplicates(subset='transcript_id')
    print(syntelogs)
    return syntelogs


def merge_blast_with_syntelogs(blast_df: pd.DataFrame, syntelogs: pd.DataFrame) -> pd.DataFrame:
    """
    Merge BLAST results with syntelog information to identify relationships.
    
    Args:
        blast_df: DataFrame with BLAST results
        syntelogs: DataFrame with syntelog information
        
    Returns:
        Merged DataFrame with BLAST and syntelog information
    """
    # Merge BLAST results with syntelogs
    merged = pd.merge(blast_df, syntelogs, left_on='query', right_on='transcript_id', how='inner')
    merged = pd.merge(merged, syntelogs, left_on='subject', right_on='transcript_id', how='inner', 
                      suffixes=('_x', '_y'))
    
    # Keep only pairs from the same syntenic group
    merged = merged[merged['Synt_id_x'] == merged['Synt_id_y']]
    
    # Create sorted haplotype combinations
    merged['haplotype_comb'] = merged[['haplotype_x', 'haplotype_y']].apply(
        lambda x: '-'.join(sorted(x)), axis=1
    )
    
    # Select relevant columns and remove duplicates
    merged = merged[[
        'haplotype_comb', 'Synt_id_x', 'identity', 'mismatch', 
        'CDS_haplotype_with_longest_annotation_x', 'CDS_length_category_x', 'query'
    ]]
    
    # Sort by identity and mismatch, then remove duplicates
    merged = (merged.sort_values(by=['identity', 'mismatch'])
              .drop_duplicates(['haplotype_comb', 'Synt_id_x'], keep='first'))
    
    return merged


def categorize_sequence_similarity(df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Categorize sequence similarity based on identity and mismatch patterns.
    
    Args:
        df: DataFrame with merged BLAST and syntelog data
        
    Returns:
        Tuple containing:
            1. DataFrame with sequence similarity categories
            2. DataFrame with identity match statistics
    """
    # Group by syntelogs and sequence similarity metrics
    grouped = (df.groupby(['Synt_id_x', 'identity', 'mismatch'])
               .size()
               .reset_index(name='counts_of_identical_combinations'))
    
    # Categorize as identical or SNP-containing
    grouped['mismatch_category'] = np.where(
        (grouped['identity'] == 100.000) & (grouped['mismatch'] == 0),
        'identical',
        'SNPs'
    )
    
    # Add rows for Synt_ids with no identical matches
    new_rows = []
    for synt_id in grouped['Synt_id_x'].unique():
        if 'identical' not in grouped[grouped['Synt_id_x'] == synt_id]['mismatch_category'].values:
            new_rows.append({
                'Synt_id_x': synt_id,
                'identity': 100.000,
                'mismatch': 0,
                'counts_of_identical_combinations': 0,
                'mismatch_category': 'identical'
            })
    
    if new_rows:
        grouped = pd.concat([grouped, pd.DataFrame(new_rows)], ignore_index=True)
    
    # Extract identical matches
    identical_matches = grouped[grouped['mismatch_category'] == 'identical'].sort_values(by='Synt_id_x')
    
    return grouped, identical_matches


def assign_allele_categories(df: pd.DataFrame) -> pd.DataFrame:
    """
    Assign allele categories based on the number of identical combinations.
    
    Args:
        df: DataFrame with identical match counts
        
    Returns:
        DataFrame with allele categories
    """
    # Map counts to meaningful categories
    category_map = {
        0: 'all_alleles_different',
        1: '2_alleles_identical',
        2: '2x2_alleles_identical',
        3: '3_alleles_identical',
        6: '4_alleles_identical'
    }
    
    df['category'] = 'no_cat'
    for count, category in category_map.items():
        df.loc[df['counts_of_identical_combinations'] == count, 'category'] = category
    
    # Remove rows with undefined categories
    df = df[df['category'] != 'no_cat']
    
    return df


def create_visualizations(data: Dict[str, pd.DataFrame], output_prefix: str) -> None:
    """
    Create and save visualizations based on the analysis results.
    
    Args:
        data: Dictionary containing analysis DataFrames
        output_prefix: Prefix for output files
    """
    # Plot allele category histogram
    identity_stats = data['identity_stats'].copy()
    identity_stats.plot(kind='bar', x='category', y='counts_per_synt_id')
    plt.title('Distribution of Allelic Categories')
    plt.xlabel('Allelic Category')
    plt.ylabel('Number of Syntenic Groups')
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_allele_categories.png', dpi=300)
    plt.close()
    
    # Plot SNP histogram
    snps_df = data['merged_with_categories'][data['merged_with_categories']['mismatch_category'] == 'SNPs']

    print(snps_df)
    
    if not snps_df.empty:
        snps_df['mismatch'] = snps_df['mismatch'].astype(int)
        plt.figure(figsize=(10, 6))
        plt.hist(snps_df['mismatch'], bins=50, alpha=0.75)
        plt.xlabel('Number of SNPs')
        plt.ylabel('Frequency')
        plt.title('Distribution of SNPs in Syntenic Genes')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(f'{output_prefix}_SNP_distribution.png', dpi=300)
        plt.close()


def main():
    """Main function to execute the analysis workflow."""
    # Parse command line arguments
    args = parse_arguments()
    
    print("=== Syntenic Gene Analysis ===")
    print(f"BLAST file: {args.blast}")
    print(f"Syntelogs file: {args.syntenic_genes}")
    print(f"Output prefix: {args.output_prefix}")
    print("-----------------------------")
    
    # Step 1: Parse BLAST results
    print("Parsing BLAST results...")
    blast_df = parse_blast_fmt6(args.blast)
    print(f"Found {len(blast_df)} BLAST alignments after filtering self-matches")
    
    # Step 2: Load syntelogs data
    print("Loading syntelogs data...")
    syntelogs = load_syntelogs(args.syntenic_genes)
    filtered_syntelogs = filter_syntelogs(syntelogs)
    print(f"Found {len(syntelogs)} syntelog entries with equal CDS lengths")
    
    # Step 3: Merge BLAST results with syntelogs
    print("Merging BLAST results with syntelogs...")

    merged_data = merge_blast_with_syntelogs(blast_df, filtered_syntelogs)
    print(f"Number of unique syntenic IDs: {len(merged_data['Synt_id_x'].unique())}")
    
    # Step 4: Categorize sequence similarity
    print("Categorizing sequence similarity...")
    similarity_data, identical_matches = categorize_sequence_similarity(merged_data)
    
    # Step 5: Assign allele categories
    print("Assigning allele categories...")
    categorized_identical = assign_allele_categories(identical_matches)
    
    # Step 6: Add categories to the main dataset
    merged_with_categories = pd.merge(
        syntelogs, 
        similarity_data[['Synt_id_x', 'mismatch_category', 'mismatch']],
        left_on ='Synt_id',
        right_on='Synt_id_x',
        how='left'
    )
    print(merged_with_categories)
    print(merged_with_categories.columns)
    # drop duplicated query
    merged_with_categories = merged_with_categories.drop_duplicates(subset='transcript_id')
    
    # Select and output relevant columns
    output_data = merged_with_categories[[
        'transcript_id', 'Synt_id',   'synteny_category',
        'CDS_haplotype_with_longest_annotation', 'CDS_length_category', 'mismatch_category', 'mismatch'
    ]]
    
    output_path = f'{args.output_prefix}_syntelog_blast_analysis.tsv'
    print(f"Writing full results to {output_path}")
    output_data.to_csv(output_path, sep='\t', index=False)
    
    # Prepare data for visualizations
    identity_group_stats = (categorized_identical.groupby(['counts_of_identical_combinations', 'category'])
                           .size()
                           .reset_index(name='counts_per_synt_id'))
    print(categorized_identical)
    # Step 7: Create visualizations
    print("Creating visualizations...")
    create_visualizations({
        'similarity_data': similarity_data,
        'merged_with_categories': merged_with_categories,
        'identity_stats': identity_group_stats
    }, args.output_prefix)
    
    print("Analysis complete!")


if __name__ == '__main__':
    main()