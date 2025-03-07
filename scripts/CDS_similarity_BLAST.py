import pandas as pd
import argparse
import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List
from pathlib import Path

def parse_blast_fmt6(blast_file: str) -> pd.DataFrame:
    """
    Parse BLAST output file in format 6 into a pandas DataFrame.
    """
    columns = ['query', 'subject', 'identity', 'length', 'mismatch', 'gap',
               'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bitscore']
   
    return pd.read_csv(blast_file, sep='\t', header=None, names=columns)

def filter_blast(blast_df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter BLAST results to remove self-matches.

    """
    return blast_df[blast_df['query'] != blast_df['subject']]

def process_syntelogs(syntelogs_file: str) -> pd.DataFrame:
    """
    Process syntelogs file and prepare for analysis.
    
    Args:
        syntelogs_file: Path to syntelogs TSV file
        
    Returns:
        Processed syntelogs DataFrame
    """
    syntelogs = pd.read_csv(
        syntelogs_file,
        sep='\t',
        usecols=['Synt_id', 'haplotype_id', 'CDS_haplotype_with_longest_annotation']
    )
    
    return syntelogs

def filter_same_length_syntelogs(syntelogs: pd.DataFrame) -> pd.DataFrame:
    """
    Filter syntelogs for equal lengths and extract haplotype information.
    
    Args:
        syntelogs: DataFrame containing syntelogs data
    """
    # Filter for equal lengths and extract haplotype information
    syntelogs = syntelogs[syntelogs['CDS_haplotype_with_longest_annotation'] == 'equal_lengths']
    syntelogs['haplotype'] = syntelogs['haplotype_id'].str.extract(r'(H\d+)')

    return syntelogs

def create_mismatch_categories(df: pd.DataFrame) -> pd.DataFrame:
    """
    Create mismatch categories based on identity and mismatch values.
    
    Args:
        df: DataFrame containing identity and mismatch columns
        
    Returns:
        DataFrame with added mismatch categories
    """
    df['mismatch_category'] = np.where(
        (df['identity'] == 100.000) & (df['mismatch'] == 0),
        'identical',
        'SNPs'
    )
    return df

def add_missing_identical_rows(df: pd.DataFrame) -> pd.DataFrame:
    """
    Some Synt_ids dont have identical matches, so  add a row with 0 mismatch and 100% identity
    and count 0.
    """
    new_rows = []
    for synt_id in df['Synt_id_x'].unique():
        if 'identical' not in df[df['Synt_id_x'] == synt_id]['mismatch_category'].values:
            new_rows.append({
                'Synt_id_x': synt_id,
                'identity': 100.000,
                'mismatch': 0,
                'counts_of_identical_combinations': 0,
                'mismatch_category': 'identical'
            })
    
    if new_rows:
        return pd.concat([df, pd.DataFrame(new_rows)], ignore_index=True)
    return df

def categorize_counts(df: pd.DataFrame) -> pd.DataFrame:
    """
    Based on the frequency of identical matche combination per Synt_id,
    add the corresponding category, as this is not intuitive from the counts.
    
    Args:
        df: DataFrame to categorize
        
    Returns:
        DataFrame with added categories
    """
    category_map = {
        0: 'all_alleles_different',
        1: '2_alleles_identical',
        2: '2x2_alleles_identical',
        3: '3_alleles_identical',
        6: '4_alleles_identical'
    }
    
    df['category'] = 'no_cat'
    for count, category in category_map.items():
        df['category'] = np.where(
            df['counts_of_identical_combinations'] == count,
            category,
            df['category']
        )
    # remove rows with no category as this is some BLAST error
    df = df[df['category'] != 'no_cat']
    return df

def plot_haplotype_histogram(df: pd.DataFrame, output_path: str) -> None:
    """
    Create and save haplotype histogram plot.
    
    Args:
        df: DataFrame containing haplotype data
        output_path: Path to save the plot
    """
    df.plot(kind='bar', x='category', y='counts_per_synt_id')
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def plot_snp_histogram(df: pd.DataFrame, output_path: str) -> None:
    """
    Create and save SNP histogram plot.
    
    Args:
        df: DataFrame containing SNP data
        output_path: Path to save the plot
    """
    df['mismatch'] = df['mismatch'].astype(int)
    counts_mismatches = df['mismatch'].value_counts()
    counts_mismatches.plot(kind='hist', bins=100)
    plt.xlabel('Number of SNPs')
    plt.ylabel('Counts')
    plt.title('Number of SNPs per Syntenic ID')
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Analyze BLAST results and syntenic genes')
    parser.add_argument('-b', '--blast', help='The blast all-vs-all', required=True)
    parser.add_argument('-s', '--synthenic_genes', help='Tsv file with syntenic genes', required=True)
    parser.add_argument('-o', '--output_prefix', help='Prefix for output files', required=True)
    
    args = parser.parse_args()
    
    
    
    # Process BLAST results
    blast_df = parse_blast_fmt6(args.blast)
    filtered_blast = filter_blast(blast_df)
    
    # Process syntelogs
    syntelogs_all = process_syntelogs(args.synthenic_genes)
    syntelogs_same_length = filter_same_length_syntelogs(syntelogs_all)
    
    # Merge data
    merged = pd.merge(filtered_blast, syntelogs, left_on='query', right_on='haplotype_id', how='inner')
    merged2 = pd.merge(merged, syntelogs, left_on='subject', right_on='haplotype_id', how='inner')
    merged2 = merged2[merged2['Synt_id_x'] == merged2['Synt_id_y']]
    
    print(f"Number of unique syntenic IDs: {len(merged2['Synt_id_x'].unique())}")
    
    # Create haplotype combinations
    merged2['haplotype_comb'] = merged2[['haplotype_x', 'haplotype_y']].apply(
        lambda x: '-'.join(sorted(x)), axis=1
    )
    

    # Process for analysis
    merged2 = merged2[['haplotype_comb', 'Synt_id_x', 'identity', 'mismatch', 'CDS_haplotype_with_longest_annotation_x', 'query']]
    merged2 = (merged2.sort_values(by=['identity', 'mismatch'])
               .drop_duplicates(['haplotype_comb', 'Synt_id_x'], keep='first'))
    
    # Group and categorize
    merged2_group = (merged2.groupby(['Synt_id_x', 'identity', 'mismatch'])
               .size()
               .reset_index(name='counts_of_identical_combinations'))
    
    merged2_group = create_mismatch_categories(merged2_group)
    merged2_group = add_missing_identical_rows(merged2_group)
    merged2_group = merged2_group.sort_values(by='Synt_id_x')
    print(merged2_group)
    # Process identical matches
    merged2_ident = merged2_group[merged2_group['mismatch_category'] == 'identical']
    merged2_ident2 = categorize_counts(merged2_ident)
    # add the category to the exploded table
    merged2 = pd.merge(merged2, merged2_ident2[['Synt_id_x', 'category']], on='Synt_id_x', how='left')
    # only keep query mismatch Synt_id_x, CDS_haplotype_with_longest_annotation_x, category
    merged2 = merged2[['query','mismatch','Synt_id_x', 'CDS_haplotype_with_longest_annotation_x', 'category']]
    # Save full results
    output_path = f'{args.output_prefix}_De_v1.functional_annotations_nodigits_genespace_syntelogs_all_exploded_blast.tsv'
    merged2.to_csv(output_path, sep='\t', index=False)
    
    merged2_ident = (merged2_ident.groupby(['counts_of_identical_combinations', 'mismatch'])
                    .size()
                    .reset_index(name='counts_per_synt_id'))
    
    merged2_ident = categorize_counts(merged2_ident)
    
    # Create plots
    plot_haplotype_histogram(merged2_ident,f'{args.output_prefix}_haplotype_histogram.png')
    
    # Process and plot SNPs
    merged2_SNPs = merged2[merged2['mismatch_category'] == 'SNPs']
    plot_snp_histogram(merged2_SNPs, f'{args.output_prefix}_SNP_histogram.png')

if __name__ == '__main__':
    main()



    # python /scratch/nadjafn/potato-allelic-orthogroups/scripts/CDS_similarity_BLAST.py -b results/06_BLAST/De_v1.functional_annotations_nodigits.txt -s results/03_GENESPACE/De_v1.functional_annotations_nodigits_genespace_syntelogs_all_exploded.tsv  -o test