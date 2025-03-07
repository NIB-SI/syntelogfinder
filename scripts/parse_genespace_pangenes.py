# load packages
import pandas as pd
import gffutils
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import argparse
import os
from typing import Dict, List, Tuple, Union, Callable


# Issue: Some transcripts get multiple synt categories then in the final output they will be duplicated
# Maybe drop these duplicates in the beginning?

def parse_pangenes(pangenes_file: str) -> pd.DataFrame:
    """Parse the tab-separated pangenes file.
    
    Args:
        pangenes_file: Path to the pangenes file
        
    Returns:
        DataFrame containing the pangenes data
    """
    return pd.read_csv(pangenes_file, sep="\t", header=0, index_col=False)


def parse_gff(gff_file: str) -> pd.DataFrame:
    """Parse the GFF file to get all transcripts.
    
    Args:
        gff_file: Path to the GFF file
        
    Returns:
        DataFrame containing the parsed GFF data with transcript IDs
    """
    gff = pd.read_csv(gff_file, sep="\t", header=None, index_col=False, comment='#')
    gff.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    gff['transcript_id'] = gff['attributes'].str.extract(r'ID=(.*?);')
    return gff


def percent_func(pct: float, allvalues: np.ndarray) -> str:
    """Format percentage labels for pie charts.
    
    Args:
        pct: The percentage value
        allvalues: Array of all values in the chart
        
    Returns:
        Formatted string with percentage and absolute count
    """
    absolute = int(pct / 100. * np.sum(allvalues))
    return "{:.1f}%\n({:d})".format(pct, absolute)


def count_comma_separated_values(df: pd.DataFrame, column_name: str) -> pd.Series:
    """Count the number of comma-separated values in a column.
    
    Args:
        df: DataFrame containing the column
        column_name: Name of the column to process
        
    Returns:
        Series with counts of comma-separated values for each row
    """
    def count_values(cell):
        if pd.isna(cell):
            return 0
        return len(str(cell).split(','))

    return df[column_name].apply(count_values)


def get_attribute_lengths(gff_file: str, attribute: str) -> pd.DataFrame:
    """Extract feature lengths from GFF file for a specific attribute type.
    
    Args:
        gff_file: Path to the GFF file
        attribute: Feature type to extract (e.g., 'exon', 'CDS')
        
    Returns:
        DataFrame with transcript IDs and their corresponding feature lengths
    """
    db_filename = 'gff.db'
    
    # Check if the gff.db file exists
    if not os.path.exists(db_filename):
        # Create a database from the GFF file
        db = gffutils.create_db(gff_file, dbfn=db_filename, keep_order=False,
                             merge_strategy='merge', sort_attribute_values=False)
    
    # Load the existing gff.db file
    db = gffutils.FeatureDB(db_filename, keep_order=True)

    # Initialize a dictionary to store transcript lengths and start coordinates
    transcript_info = {}

    # Iterate over all features of the specified type in the GFF file
    for feature in db.features_of_type(attribute):
        parent_id = feature.attributes['Parent'][0]
        length = feature.end - feature.start + 1
        
        if parent_id not in transcript_info:
            # Get the parent transcript feature
            parent = db[parent_id]
            transcript_info[parent_id] = {
                f'{attribute}_ref_length': 0,
                f'{attribute}_parent_start': parent.start
            }
        
        transcript_info[parent_id][f'{attribute}_ref_length'] += length
    
    # Convert dictionary to DataFrame
    df = pd.DataFrame.from_dict(transcript_info, orient='index')
    df['transcript_id'] = df.index
    return df


def make_barplot(df: pd.DataFrame, attribute: str, ax: plt.Axes, haplotype_labels: List[str]) -> None:
    """Create a barplot for length categories on the given axes.
    
    Args:
        df: DataFrame containing the data
        attribute: Feature type to plot
        ax: Matplotlib axes to plot on
        haplotype_labels: List of haplotype labels
    """
    # Define custom order for categories
    custom_order = [
        'less_1%_difference', 
        'more_1%_difference', 
        'more_5%_difference', 
        'more_10%_difference', 
        'more_20%_difference'
    ]
    
    # Convert to categorical for consistent ordering
    df[f'{attribute}_length_category'] = pd.Categorical(
        df[f'{attribute}_length_category'], 
        categories=custom_order, 
        ordered=True
    )
    
    # Create list of possible haplotype labels plus 'equal_lengths'
    haplotype_categories = haplotype_labels + ['equal_lengths']
    
    # Change the haplotype_with_longest_annotation to a categorical variable
    df[f'{attribute}_haplotype_with_longest_annotation'] = pd.Categorical(
        df[f'{attribute}_haplotype_with_longest_annotation'], 
        categories=haplotype_categories, 
        ordered=True
    )
    
    # Create the barplot on the given axes
    sns.histplot(
        data=df, 
        x=f'{attribute}_length_category', 
        hue=f'{attribute}_haplotype_with_longest_annotation', 
        multiple='stack', 
        ax=ax
    )
    
    ax.set_xlabel(f'{attribute} Length Category')
    ax.set_ylabel('Count')
    ax.set_ylim(0, 15000)
    ax.tick_params(axis='x', rotation=90)  # Rotate x-axis labels


def generate_combined_plots(df_list: List[pd.DataFrame], attributes: List[str], output_prefix: str) -> None:
    """Generate combined barplots for multiple attributes.
    
    Args:
        df_list: List of DataFrames containing data for each attribute
        attributes: List of feature types to plot
        output_prefix: Prefix for the output file
    """
    # Set up the subplots (1 row, n columns for each attribute)
    fig, axes = plt.subplots(1, len(attributes), figsize=(15, 5), sharey=True)

    # If only one plot, axes won't be a list, so make it iterable
    if len(attributes) == 1:
        axes = [axes]

    # Loop through attributes and create a plot for each
    for i, (df, attribute) in enumerate(zip(df_list, attributes)):
        make_barplot(df, attribute, axes[i], get_haplotype_labels(df))
        if i > 0:
            axes[i].legend_.remove()  # Remove individual legends

    # Adjust layout to prevent overlap and save the combined plot
    plt.tight_layout()  # Give space for the legend
    plt.savefig(f'{output_prefix}_combined_barplots.png')
    plt.close(fig)  # Close figure to free memory


def get_haplotype_columns(df: pd.DataFrame) -> List[str]:
    """Get a list of haplotype columns from the DataFrame.
    
    Args:
        df: DataFrame containing haplotype columns
    
    Returns:
        List of haplotype column names
    """
    return [col for col in df.columns if col.startswith('hap')]


def get_haplotype_labels(df: pd.DataFrame) -> List[str]:
    """Get haplotype labels based on column names.
    
    Args:
        df: DataFrame with attribute columns containing haplotype labels
    
    Returns:
        List of haplotype labels (e.g., ['1G', '2G', '3G', '4G'])
    """
    # Extract column names that contain 'length_' followed by a number and 'G'
    columns = [col for col in df.columns if '_length_' in col and col.endswith('G')]
    
    # Extract the haplotype labels (e.g., '1G', '2G')
    labels = [col.split('_')[-1] for col in columns]
    
    return sorted(labels)


def add_synteny_category(pangenes: pd.DataFrame) -> pd.DataFrame:
    """Add synteny categories to the pangenes DataFrame.
    
    Args:
        pangenes: DataFrame containing the pangenes data
        
    Returns:
        DataFrame with added synteny categories and pivoted structure
    """
    # Get haplotype columns
    haplotype_cols = get_haplotype_columns(pangenes)
    ploidy = len(haplotype_cols)
    
    # Count genes in each haplotype
    for hap in haplotype_cols:
        pangenes[f'{hap}_count'] = count_comma_separated_values(pangenes, hap).astype(str) + hap
    
    # Add synteny classification
    pangenes['true_synteny'] = 'synteny'
    
    # Check for special characters in any haplotype column
    has_special_chars = False
    for hap in haplotype_cols:
        has_special_chars = has_special_chars | pangenes[hap].str.contains('\*') | pangenes[hap].str.contains('\+')
    
    pangenes.loc[has_special_chars, 'true_synteny'] = 'no_synteny'

    # Add a synteny ID column to group transcripts in the same synteny group
    pangenes['Synt_id'] = 'Synt_id_' + pangenes.index.astype(str)
    
    # Create synteny category combining all haplotype information
    synteny_category_parts = [pangenes[f'{hap}_count'] for hap in haplotype_cols]
    synteny_category_parts = pd.DataFrame(synteny_category_parts).T  # Transpose to get correct column format
    pangenes['synteny_category'] = synteny_category_parts.astype(str).agg('_'.join, axis=1) + '_' + pangenes['true_synteny']
    # Create clean copy without special characters
    pangenes_gene = pangenes.map(lambda x: x.replace("+", "").replace("*", "") if isinstance(x, str) else x)
    
    # Combine all haplotype columns
    syntenic_genes_parts = [pangenes_gene[hap] for hap in haplotype_cols]
    syntenic_genes_parts = pd.DataFrame(syntenic_genes_parts).T
    pangenes_gene['syntenic_genes'] = syntenic_genes_parts.astype(str).agg(','.join, axis=1)
    
    # Pivot the table for easier analysis
    pangenes_pivot = pd.melt(
        pangenes_gene, 
        id_vars=['Synt_id', 'synteny_category', 'syntenic_genes'], 
        value_vars=haplotype_cols, 
        var_name='haplotype', 
        value_name='transcript_id'
    )
    
    # Drop rows where transcript_id is NaN
    pangenes_pivot = pangenes_pivot.dropna(subset='transcript_id')

    # drop duplicated rows
    pangenes_pivot = pangenes_pivot.drop_duplicates(['transcript_id', 'synteny_category'])

    
    # Split comma-separated transcript IDs and expand the DataFrame
    pangenes_pivot['transcript_id'] = pangenes_pivot['transcript_id'].str.split(',')
    pangenes_pivot = pangenes_pivot.explode('transcript_id')

    return pangenes_pivot


def merge_pangenes_gff(pangenes_pivot: pd.DataFrame, gff: pd.DataFrame) -> pd.DataFrame:
    """Merge pangenes and GFF DataFrames.
    
    Args:
        pangenes_pivot: Pivoted pangenes DataFrame
        gff: GFF DataFrame
        
    Returns:
        Merged DataFrame
    """
    # Select only mRNA from gff
    gff_mrna = gff[gff['type'] == 'mRNA'].copy()
    
    # Ensure transcript_id columns have the same type
    gff_mrna['transcript_id'] = gff_mrna['transcript_id'].astype(str)
    pangenes_pivot['transcript_id'] = pangenes_pivot['transcript_id'].astype(str)
    
    # Merge the DataFrames
    gff_pangenes = pd.merge(gff_mrna, pangenes_pivot, on='transcript_id', how='left')
    
    # Convert synteny_category to string
    gff_pangenes['synteny_category'] = gff_pangenes['synteny_category'].astype(str)
    
    return gff_pangenes


def make_pie_chart(gff_pangenes: pd.DataFrame, output_prefix: str) -> None:
    """Create a pie chart showing distribution of synteny categories.
    
    Args:
        gff_pangenes: Merged GFF and pangenes DataFrame
        output_prefix: Prefix for the output file
    """
    # Select only mRNA rows
    mrna_data = gff_pangenes[gff_pangenes['type'] == 'mRNA'].copy()
    
    # Sort by synteny category
    mrna_data = mrna_data.sort_values('synteny_category')
    
    # Count genes in each synteny category
    synt_counts = mrna_data['synteny_category'].value_counts()
    
    # Get top 5 categories and group the rest as "other"
    synt_counts_top = synt_counts[:5].copy()
    synt_counts_top['other'] = synt_counts[5:].sum()
    
    # Sort alphabetically
    synt_counts_top.sort_index(inplace=True)
    
    # Set up explode values for pie chart
    explode = tuple([0.1 * (6-i) for i in range(len(synt_counts_top))])
    
    # Create figure
    plt.figure(figsize=(10, 8))
    
    # Create pie chart
    synt_counts_top.plot.pie(
        startangle=90, 
        explode=explode, 
        autopct=lambda pct: percent_func(pct, synt_counts)
    )
    
    # Save the chart
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_pie_chart.png')
    plt.close()  # Close figure to free memory


def check_length_values(row: pd.Series, percent: float) -> bool:
    """Check if the difference between min and max values exceeds the specified percentage.
    
    Args:
        row: Series containing length values
        percent: Percentage threshold
        
    Returns:
        True if difference exceeds threshold, False otherwise
    """
    min_value = row.min()
    max_value = row.max()

    # Calculate the percentage cutoff
    cutoff = min_value * percent * 0.01

    # Check if max value exceeds min value plus cutoff
    return max_value >= min_value + cutoff


def add_length_category(df: pd.DataFrame, attribute: str) -> pd.DataFrame:
    """Categorize length differences for transcripts.
    
    Args:
        df: DataFrame containing length data
        attribute: Feature type to categorize
        
    Returns:
        DataFrame with added length categories
    """
    # Select length columns for all haplotypes
    length_cols = [col for col in df.columns if f'{attribute}_length_' in col and col.endswith('G')]
    
    df_subset = df[length_cols].astype(int)
    
    # Initialize category column
    df[f'{attribute}_length_category'] = 'unclassified'
    
    # Apply categorization rules
    df.loc[~df_subset.apply(lambda x: check_length_values(x, 1), axis=1), 
           f'{attribute}_length_category'] = 'less_1%_difference'
    
    df.loc[df_subset.apply(lambda x: check_length_values(x, 1), axis=1) & 
           ~df_subset.apply(lambda x: check_length_values(x, 5), axis=1), 
           f'{attribute}_length_category'] = 'more_1%_difference'
    
    df.loc[df_subset.apply(lambda x: check_length_values(x, 5), axis=1) & 
           ~df_subset.apply(lambda x: check_length_values(x, 10), axis=1), 
           f'{attribute}_length_category'] = 'more_5%_difference'
    
    df.loc[df_subset.apply(lambda x: check_length_values(x, 10), axis=1) & 
           ~df_subset.apply(lambda x: check_length_values(x, 20), axis=1), 
           f'{attribute}_length_category'] = 'more_10%_difference'
    
    df.loc[df_subset.apply(lambda x: check_length_values(x, 20), axis=1), 
           f'{attribute}_length_category'] = 'more_20%_difference'
    
    return df


def add_length_differences_percent(df: pd.DataFrame, attribute: str) -> pd.DataFrame:
    """Calculate percentage differences in lengths between haplotypes.
    
    Args:
        df: DataFrame containing length data
        attribute: Feature type to analyze
        
    Returns:
        DataFrame with added difference metrics
    """
    # Select length columns
    length_cols = [col for col in df.columns if f'{attribute}_length_' in col and col.endswith('G')]
    
    df_subset = df[length_cols].astype(int)
    
    # Calculate max difference
    df[f'{attribute}_max_difference'] = df_subset.max(axis=1) - df_subset.min(axis=1)
    
    # Calculate percentage difference
    df[f'{attribute}_percent_difference'] = (df[f'{attribute}_max_difference'] / 
                                           df_subset.min(axis=1) * 100)
    
    return df


def add_longest_transcript(df: pd.DataFrame, attribute: str) -> pd.DataFrame:
    """Identify the haplotype with the longest transcript for each synteny group.
    
    Args:
        df: DataFrame containing length data
        attribute: Feature type to analyze
        
    Returns:
        DataFrame with added information about longest transcript
    """
    # Select length columns
    length_cols = [col for col in df.columns if f'{attribute}_length_' in col and col.endswith('G')]
    
    # Find column with maximum value for each row
    df[f'{attribute}_haplotype_with_longest_annotation'] = df[length_cols].idxmax(axis=1)
    
    # Check if all lengths are equal
    mask = (df[length_cols].nunique(axis=1) == 1)
    df.loc[mask, f'{attribute}_haplotype_with_longest_annotation'] = 'equal_lengths'
    
    # Clean up column names
    df[f'{attribute}_haplotype_with_longest_annotation'] = (
        df[f'{attribute}_haplotype_with_longest_annotation']
        .str.replace(f'{attribute}_length_', '')
    )
    
    return df


def filter_all_syntelogs(pangenes: pd.DataFrame, target_category: str = None) -> pd.DataFrame:
    """Filter pangenes to include only specified synteny category.
    
    Args:
        pangenes: DataFrame containing pangenes data
        target_category: Synteny category to filter for. If None, will construct a
                        filter for syntelogs with 1 gene per haplotype.
                        
    Returns:
        Filtered DataFrame
    """
    if target_category is None:
        # Get haplotype columns
        haplotype_cols = [col for col in pangenes.columns if col.startswith('hap')]
        ploidy = len(haplotype_cols)
        
        # Construct filter category for syntelogs with 1 gene per haplotype
        filter_parts = [f"1{hap}" for hap in haplotype_cols]
        target_category = '_'.join(filter_parts) + '_synteny'
    
    # Select only rows with the specified category
    filtered = pangenes[pangenes['synteny_category'] == target_category].copy()
    return filtered


def add_length_to_syntelogs(pangenes: pd.DataFrame, ref_lengths: pd.DataFrame, attribute: str) -> pd.DataFrame:
    """Add length information to syntelog groups.
    
    Args:
        pangenes: DataFrame containing pangenes data
        ref_lengths: DataFrame with reference lengths
        attribute: Feature type to analyze
        
    Returns:
        DataFrame with length information for syntelog groups
    """
    # Merge the two dataframes
    df_synt_lengths = pd.merge(pangenes, ref_lengths, on='transcript_id', how='inner')

    # Group transcript IDs by synteny ID
    transcript_ids = (df_synt_lengths.groupby('Synt_id')['transcript_id']
                    .apply(list).reset_index())

    # Get unique haplotypes
    haplotypes = df_synt_lengths['haplotype'].unique()
    
    # Pivot data to have one row per synteny group
    df_synt_pivot = df_synt_lengths.pivot(
        index='Synt_id', 
        columns=['haplotype'], 
        values=[f'{attribute}_ref_length', f'{attribute}_parent_start']
    )
    
    # Flatten the columns
    df_synt_pivot.columns = ['_'.join(col) for col in df_synt_pivot.columns]
    
    # Merge with transcript IDs
    df_synt_pivot = pd.merge(df_synt_pivot, transcript_ids, on='Synt_id', how='inner')

    # Rename columns for clarity based on haplotype naming convention (e.g., hap1 -> 1G)
    renamed_columns = ['Synt_id']
    
    for hap in sorted(haplotypes):
        hap_num = hap.replace('hap', '')
        renamed_columns.append(f'{attribute}_length_{hap_num}G')
    
    for hap in sorted(haplotypes):
        hap_num = hap.replace('hap', '')
        renamed_columns.append(f'{attribute}_parent_start_{hap_num}G')
    
    renamed_columns.append('transcript_id')
    
    # Rename columns
    df_synt_pivot.columns = renamed_columns
    
    # Drop rows with missing values
    df_synt_pivot = df_synt_pivot.dropna()
    
    # Add length categories and differences
    df_length = add_length_category(df_synt_pivot, attribute)
    df_length = add_length_differences_percent(df_length, attribute)
    df_length_cat = add_longest_transcript(df_length, attribute)

    return df_length_cat


def main():
    """Main function to execute the analysis workflow."""
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Parse the pan-genome file and GFF file')

    parser.add_argument('-p', '--pangenes', help='The pan-genome file', required=True)
    parser.add_argument('-g', '--gff', help='The GFF file', required=True)
    parser.add_argument('-o', '--output', help='Output prefix', required=True)
    parser.add_argument('-s', '--syntelogs_category', help='Target syntelog category to filter for (default: all haplotypes with 1 gene)', default=None)

    args = parser.parse_args()

    # Parse input files
    print("Parsing pan-genome file...")
    pangenes = parse_pangenes(args.pangenes)
    
    print("Parsing GFF file...")
    gff = parse_gff(args.gff)

    # Determine ploidy level
    haplotype_cols = get_haplotype_columns(pangenes)
    ploidy = len(haplotype_cols)
    print(f"Detected ploidy level: {ploidy}x (haplotype columns: {','.join(haplotype_cols)})")

    # Generate synteny categories
    print("Adding synteny categories...")
    pangenes_pivot = add_synteny_category(pangenes)
    
    # Merge pangenes and GFF data
    print("Merging pangenes and GFF data...")
    gff_pangenes = merge_pangenes_gff(pangenes_pivot, gff)
    
    # Create pie chart
    print("Creating pie chart...")
    make_pie_chart(gff_pangenes, args.output)

    # Filter for syntelogs based on category
    print("Filtering syntelogs...")
    syntelogs = filter_all_syntelogs(pangenes_pivot, args.syntelogs_category)

    # Process each attribute type
    print("Processing feature lengths...")
    attributes = ['exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR']
    syntelogs_lengths_list = []
    
    for attribute in attributes:
        print(f"Processing {attribute} features...")
        # Get feature lengths from GFF
        ref_lengths = get_attribute_lengths(args.gff, attribute)
        
        # Add length information to syntelogs
        syntelogs_lengths = add_length_to_syntelogs(syntelogs, ref_lengths, attribute)
        
        syntelogs_lengths_list.append(syntelogs_lengths)
    
    # Generate combined plots
    print("Generating combined plots...")
    generate_combined_plots(syntelogs_lengths_list, attributes, args.output)
    # explode
    syntelogs_lengths_cds = syntelogs_lengths_list[1].explode('transcript_id')
    # drop duplicated rows
    syntelogs_lengths_cds = syntelogs_lengths_cds.drop_duplicates('transcript_id')
    # Merge transcript information with length categories
    print("Merging and saving final output...")
    print(syntelogs_lengths_cds['transcript_id'])
    if len(syntelogs_lengths_list) > 1:  # Check if we have CDS info
        pangenes_pivot = pd.merge(
            pangenes_pivot, 
            syntelogs_lengths_cds,  # CDS info is at index 1
            on='transcript_id', 
            how='left'
        )

        # Set transcript ID as index
        pangenes_pivot.set_index('transcript_id', inplace=True)
        
        # Select relevant columns for output
        final_output = pangenes_pivot[[
            'Synt_id_x', 'synteny_category', 'syntenic_genes', 'haplotype', 
            'CDS_length_category', 'CDS_percent_difference', 'CDS_haplotype_with_longest_annotation'
        ]]

        final_output.rename(columns={'Synt_id_x': 'Synt_id'}, inplace=True)
        # print duplicated rows
        print(final_output[final_output.duplicated()])
        # drop duplicated rows
        final_output = final_output.drop_duplicates()
        # print rows where transcript_id is duplicated
        print(final_output[final_output.index.duplicated()])
        # Save to TSV file
        final_output.to_csv(f'{args.output}_categories.tsv', sep='\t', index=True)
    
    print("Analysis complete!")


if __name__ == '__main__':
    main()



# python /DKED/scratch/nadjafn/potato-allelic-orthogroups/scripts/parse_pangenes.py --pangenes De_v1.functional_annotations_nodigits_genespace.tsv                                               --gff De_v1.functional_annotations_nodigits.longest_isoforms.gff                                                          --output De_v1.functional_annotations_nodigits_genespace -s 1hap1_1hap2_1hap3_1hap4_synteny