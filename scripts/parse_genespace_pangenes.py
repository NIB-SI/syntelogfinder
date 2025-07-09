# load packages
import pandas as pd
import gffutils
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import argparse
import os
from typing import Dict, List, Tuple, Union, Callable, Optional, Set
import re
from collections import defaultdict, deque

# Issue: Some transcripts get multiple synt categories then in the final output they will be duplicated
# Maybe drop these duplicates in the beginning?

def parse_pangenes(pangenes_file: str) -> pd.DataFrame:
    """Parse the tab-separated pangenes file.

    Args:
        pangenes_file: Path to the pangenes file

    Returns:
        DataFrame containing the pangenes data
    """
    df = pd.read_csv(pangenes_file, sep="\t", header=0, index_col=False, na_values=['NA', 'N/A', '', 'NaN'])
    return df.replace(', ', ',', regex=True)

def identify_overlapping_rows(df: pd.DataFrame, na_values: Optional[List] = None) -> Dict[int, Set[int]]:
    """
    Identify groups of rows that have overlapping gene IDs.

    Parameters:
    -----------
    df : pd.DataFrame
        The input dataframe
    na_values : list, optional
        List of values to treat as NA besides pandas defaults

    Returns:
    --------
    Dict[int, Set[int]]
        Dictionary mapping row indices to sets of connected row indices (including itself)
    """
    if na_values is None:
        na_values = ['NA', 'N/A', '', None]


    gene_pattern = r'(\S+)$'


    # Map gene IDs to the rows they appear in
    gene_to_rows = defaultdict(set)

    # Build the gene-to-rows mapping
    for row_idx, row in df.iterrows():
        for col in df.columns:
            value = row[col]

            # Skip NA values
            if pd.isna(value) or value in na_values:
                continue

            # Handle multiple comma-separated gene IDs in a cell
            if isinstance(value, str):
                # Find all gene IDs in this cell
                matches = re.findall(gene_pattern, value)
                for gene_id in matches:
                    gene_to_rows[gene_id].add(row_idx)

    # Create a graph of connected rows
    row_connections = defaultdict(set)

    # Connect rows that share at least one gene
    for gene, rows in gene_to_rows.items():
        if len(rows) > 1:
            # All rows with this gene are connected to each other
            for row in rows:
                row_connections[row].update(rows)

    # Find connected components using BFS
    components = {}
    visited = set()

    for start_row in row_connections:
        if start_row in visited:
            continue

        # BFS to find all connected rows
        queue = deque([start_row])
        component = set()

        while queue:
            row = queue.popleft()
            if row in visited:
                continue

            visited.add(row)
            component.add(row)

            # Add all connected rows
            for connected_row in row_connections[row]:
                if connected_row not in visited:
                    queue.append(connected_row)

        # Store this connected component
        for row in component:
            components[row] = component

    # Add singleton rows (rows not connected to any other rows)
    for row_idx in df.index:
        if row_idx not in components:
            components[row_idx] = {row_idx}

    return components

def merge_overlapping_rows(df: pd.DataFrame, na_values: Optional[List] = None,
                          print_operations: bool = True) -> pd.DataFrame:
    """
    Merge rows that have overlapping gene IDs into consolidated rows.

    Parameters:
    -----------
    df : pd.DataFrame
        The input dataframe
    na_values : list, optional
        List of values to treat as NA besides pandas defaults
    print_operations : bool, default=True
        Whether to print information about the merge operations

    Returns:
    --------
    pd.DataFrame
        A new dataframe with merged rows
    """
    if na_values is None:
        na_values = ['NA', 'N/A', '', None]

    # Identify groups of rows that need to be merged
    components = identify_overlapping_rows(df, na_values)

    # Unique components (we only need one representative per component)
    unique_components = {}
    for component in components.values():
        # Use the minimum row index as the representative
        rep = min(component)
        if rep not in unique_components:
            unique_components[rep] = component

    if print_operations:
        num_groups = len(unique_components)
        num_affected_rows = sum(len(comp) for comp in unique_components.values())
        num_merged = num_affected_rows - num_groups

        print(f"Found {num_groups} groups of overlapping rows.")
        print(f"Merging {num_merged} rows into their respective groups.\n")

    # Pattern to match gene IDs including special characters
    #gene_pattern = r'(Soltu\.Des\.v1_[A-Za-z0-9\.]+\d+[\*\+]?)'
    #gene_pattern = r'[A-Za-z0-9]+(?:\.[A-Za-z0-9]+)?_[A-Za-z0-9]+(?:\.\d+[\*\+]?)'
    # More specific but still general pattern
    # gene_pattern = r'([A-Za-z0-9\.]+_C\d+H\d+G\d+\.\d+[\*\+]?)'
    #gene_pattern = r'([A-Za-z0-9]+[._][A-Za-z0-9_]+\d+[A-Za-z0-9_.]*\*?[\*\+]?)'
    gene_pattern = r'(\S+)$'
    ## here is the problem with the regex, it is not matching the gene IDs correctly

    # Create the merged dataframe
    result_data = []

    # Process each component
    for rep, component in unique_components.items():
        if len(component) == 1:
            # Single row, no merging needed
            result_data.append(df.loc[rep].to_dict())
            continue

        if print_operations:
            print(f"Merging rows: {', '.join(map(str, sorted(component)))}")

        # Create a new merged row
        merged_row = {}

        # For each column, merge all gene IDs from the component rows
        for col in df.columns:
            gene_ids = set()

            for row_idx in component:
                value = df.at[row_idx, col]

                # Skip NA values
                if pd.isna(value) or value in na_values:
                    continue

                # Extract gene IDs including special characters
                if isinstance(value, str):
                    matches = re.findall(gene_pattern, value)
                    gene_ids.update(matches)

            # Ensure no empty strings in gene_ids
            gene_ids.discard('')

            # Create the merged cell value
            if gene_ids:
                merged_row[col] = ','.join(sorted(gene_ids))  # No need for set() again
            else:
                merged_row[col] = pd.NA


        result_data.append(merged_row)

    # Create the result dataframe
    result_df = pd.DataFrame(result_data, columns=df.columns)



    if print_operations:
        print(f"\nOriginal dataframe had {len(df)} rows.")
        print(f"New dataframe has {len(result_df)} rows.")
        print(f"Reduced by {len(df) - len(result_df)} rows.")

    return result_df

def remove_duplicates(df):
    def process_cell(cell):
        if not isinstance(cell, str):  # Skip NaN values
            return cell

        gene_ids = cell.split(',')

        # Using a dictionary for quick lookups
        unique_genes = {}
        for gene in gene_ids:
            base_gene = re.sub(r'[\*\+]', '', gene)  # Remove * and + for comparison
            # Keep only the one with * or + if available
            if base_gene not in unique_genes or '*' in gene or '+' in gene:
                unique_genes[base_gene] = gene

        return ','.join(unique_genes.values())

    # Use Pandas' vectorized `map()` instead of `applymap()` for speed
    return df.map(process_cell)

def parse_gff(gff_file: str) -> pd.DataFrame:
    """Parse the GFF file to get all transcripts.
    Args:
        gff_file: Path to the GFF file
    Returns:
        DataFrame containing the parsed GFF data with transcript IDs
    """
    gff = pd.read_csv(gff_file, sep="\t", header=None, index_col=False, comment='#')
    gff.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']

    # Extract transcript_id
    gff['transcript_id'] = np.where(
        gff['attributes'].str.contains('Parent='),
        gff['attributes'].str.extract(r'Parent=([^;]+)')[0],
        gff['attributes'].str.extract(r'ID=([^;]+)')[0]
    )

    # First, extract gene_id for mRNA/transcript entries
    gff['gene_id'] = np.where(
        gff['type'].str.contains('mRNA|transcript'),
        gff['attributes'].str.extract(r'geneID=([^;]+)')[0],
        ''
    )

    # Create a mapping of transcript_id to gene_id for mRNA/transcript entries
    transcript_to_gene = gff[gff['gene_id'] != ''].set_index('transcript_id')['gene_id'].to_dict()

    # Propagate gene_id to all rows with the same transcript_id
    gff['gene_id'] = gff['transcript_id'].map(transcript_to_gene).fillna('')

    print(gff)
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
    """Extract feature lengths from a GFF file for a specific attribute type.

    Args:
        gff_file: Path to the GFF file.
        attribute: Feature type to extract (e.g., 'exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR').

    Returns:
        DataFrame with transcript IDs and their corresponding feature lengths.
    """
    db_filename = 'gff.db'

    # Create a database if it doesn't exist
    if not os.path.exists(db_filename):
        print(f"Creating GFF database: {db_filename}")
        db = gffutils.create_db(gff_file, dbfn=db_filename, keep_order=False,
                                merge_strategy="error", disable_infer_transcripts=True, verbose = True )
    else:
        print(f"Using existing GFF database: {db_filename}")

    # Load the GFF database
    db = gffutils.FeatureDB(db_filename, keep_order=True)

    transcript_info = {}

    # Collect all transcript IDs in the GFF file
    all_transcripts = {feature.id for feature in db.features_of_type("mRNA")}

    # Iterate over all features of the specified type
    for feature in db.features_of_type(attribute):
        # Handle missing Parent attribute safely
        parent_id = feature.attributes.get('Parent', [None])[0]

        if parent_id is None:
            parent_id = feature.attributes.get('ID', [None])[0]
            #continue  # Skip if no parent ID found

        length = feature.end - feature.start + 1

        if parent_id not in transcript_info:
            # Try retrieving the parent feature safely
            try:
                parent = db[parent_id]
                parent_start = parent.start
            except gffutils.FeatureNotFoundError:
                parent_start = None  # Handle missing parent case

            transcript_info[parent_id] = {
                f'{attribute}_ref_length': 0,
                f'{attribute}_parent_start': parent_start
            }

        transcript_info[parent_id][f'{attribute}_ref_length'] += length

    # Ensure all transcripts are in the final output
    for transcript_id in all_transcripts:
        if transcript_id not in transcript_info:
            transcript_info[transcript_id] = {
                f'{attribute}_ref_length': 0,
                f'{attribute}_parent_start': None  # Could also be 0 if preferred
            }

    # Convert dictionary to DataFrame
    df = pd.DataFrame.from_dict(transcript_info, orient='index')
    df.reset_index(inplace=True)
    df.rename(columns={'index': 'transcript_id'}, inplace=True)
    print(df)

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
        'not_annotated',
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
    ax.set_ylim(0, 25000)
    ax.tick_params(axis='x', rotation=90)  # Rotate x-axis labels


def generate_combined_plots(df_list: List[pd.DataFrame], attributes: List[str], output_prefix: str) -> None:
    """Generate combined barplots for multiple attributes.

    Args:
        df_list: List of DataFrames containing data for each attribute
        attributes: List of feature types to plot
        output_prefix: Prefix for the output file
    """
    # Set up the subplots (1 row, n columns for each attribute)
    fig, axes = plt.subplots(1, len(attributes), figsize=(13, 5), sharey=True)

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
    plt.savefig(f'{output_prefix}_combined_barplots.svg')
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
    print(pangenes)
    # Count genes in each haplotype
    for hap in haplotype_cols:
        pangenes[f'{hap}_count'] = count_comma_separated_values(pangenes, hap).astype(str) + hap

    # Add synteny classification
    pangenes['true_synteny'] = 's'

    # Check for special characters in any haplotype column
    has_special_chars = False
    for hap in haplotype_cols:
        has_special_chars = has_special_chars | pangenes[hap].str.contains('\*') | pangenes[hap].str.contains('\+')

    pangenes.loc[has_special_chars, 'true_synteny'] = 'no_s'

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
    pangenes_gene['syntenic_genes'] = syntenic_genes_parts.astype(str).agg(
        lambda x: ','.join(v for v in x if pd.notna(v) and v not in ['NA', 'N/A', '']), axis=1
    )
    print(pangenes_gene)
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

    # Split comma-separated transcript IDs and expand the DataFrame
    pangenes_pivot['transcript_id'] = pangenes_pivot['transcript_id'].str.split(',')
    pangenes_pivot = pangenes_pivot.explode('transcript_id')

    # Drop the rows where transcript_id is NA or empty
    pangenes_pivot = pangenes_pivot.dropna(subset=['transcript_id'])
    pangenes_pivot = pangenes_pivot[pangenes_pivot['transcript_id'].str.strip() != '']

    # Remove duplicate transcript_ids by keeping first occurrence
    # This ensures each transcript only appears once in the resulting dataframe
    pangenes_pivot = pangenes_pivot.drop_duplicates(subset=['transcript_id'])

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
    print(gff_pangenes)
    # Convert synteny_category to string
    gff_pangenes['synteny_category'] = gff_pangenes['synteny_category'].astype(str)

    return gff_pangenes


def make_pie_chart(gff_pangenes: pd.DataFrame, syntelogs_category: str, output_prefix: str) -> None:
    """Create a pie chart showing distribution of synteny categories.
    Args:
    gff_pangenes: Merged GFF and pangenes DataFrame
    output_prefix: Prefix for the output file
    """
    # Select only mRNA rows
    #mrna_data = gff_pangenes[gff_pangenes['type'] == 'transcript'].copy()
    mrna_data = gff_pangenes[gff_pangenes['type'] == 'mRNA'].copy()
    # Sort by synteny category
    mrna_data = mrna_data.sort_values('synteny_category')
    # Count genes in each synteny category
    synt_counts = mrna_data['synteny_category'].value_counts()
    print(synt_counts)

    # Get top 5 categories and group the rest as "other"
    synt_counts_top = synt_counts[:7].copy()
    synt_counts_top['other'] = synt_counts[7:].sum()
    # Sort alphabetically
    synt_counts_top.sort_index(inplace=True)

    # Set up colors: red for the specified category, gray for everything else
    colors = ['#D3D3D3'] * len(synt_counts_top)  # Default gray for all slices

    # Find the index of the syntelogs category to highlight (if it exists)
    if syntelogs_category in synt_counts_top.index:
        highlight_idx = list(synt_counts_top.index).index(syntelogs_category)
        colors[highlight_idx] = '#FF0000'  # Red for the highlighted category

    # Set up explode values for pie chart
    explode = tuple([0.1 * (7-i) for i in range(len(synt_counts_top))])
    # Create figure
    plt.figure(figsize=(7, 7))
    # Create pie chart
    synt_counts_top.plot.pie(
    startangle=90,
    explode=explode,
    autopct=lambda pct: percent_func(pct, synt_counts),
    colors = colors
    )
    # Save the chart
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_pie_chart.svg', bbox_inches='tight')
    plt.close() # Close figure to free memory


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
    # NA to 0
    df[length_cols] = df[length_cols].fillna(0)
    print(df)
    # Convert length columns to integers
    df_subset = df[length_cols].astype(int)
    print(df_subset)

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

    # Assign 'not_annotated' if all lengths are 0
    df.loc[(df_subset == 0).all(axis=1), f'{attribute}_length_category'] = 'not_annotated'
    # print the not_annotated rows
    #print(df.loc[(df_subset == 0).all(axis=1)])

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
    # df_synt_pivot = df_synt_pivot.dropna(subset=[f'{attribute}_length_1G'])

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

    print("Dropping subset rows...")
    #pangenes = drop_subset_rows(pangenes)
    print(pangenes)

    print("Merging overlapping rows...")
    pangenes = merge_overlapping_rows(pangenes, print_operations=False)
    print(pangenes)
    pangenes = remove_duplicates(pangenes)
    # save to TSV file
    pangenes.to_csv(f'{args.output}_merged_overlapping.tsv', sep='\t', index=False)
    print("Parsing GFF file...")
    gff = parse_gff(args.gff)

    # Determine ploidy level
    haplotype_cols = get_haplotype_columns(pangenes)
    ploidy = len(haplotype_cols)
    print(f"Detected ploidy level: {ploidy}x (haplotype columns: {','.join(haplotype_cols)})")

    # Generate synteny categories
    print("Adding synteny categories...")
    pangenes_pivot = add_synteny_category(pangenes)
    print(pangenes_pivot)
    # Merge pangenes and GFF data
    print("Merging pangenes and GFF data...")
    gff_pangenes = merge_pangenes_gff(pangenes_pivot, gff)

    # Create pie charts
    print("Creating pie chart...")
    make_pie_chart(gff_pangenes, args.syntelogs_category, args.output)

    # Filter for syntelogs based on category
    print("Filtering syntelogs...")
    syntelogs = filter_all_syntelogs(pangenes_pivot, args.syntelogs_category)

    # Process each attribute type
    print("Processing feature lengths...")
    # attributes = ['exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR']
    attributes = ['exon', 'CDS']
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
    # syntelogs_lengths_cds = syntelogs_lengths_cds.drop_duplicates('transcript_id')
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
        pangenes_pivot = pangenes_pivot.drop_duplicates('transcript_id')
        # Set transcript ID as index
        pangenes_pivot.set_index('transcript_id', inplace=True)

        # Select relevant columns for output
        final_output = pangenes_pivot[[
            'Synt_id_x', 'synteny_category', 'syntenic_genes', 'haplotype',
            'CDS_length_category', 'CDS_percent_difference', 'CDS_haplotype_with_longest_annotation'
        ]]
        # add the gene_id from gff_pangenes to the final output
        final_output = pd.merge(
            final_output,
            gff_pangenes[['transcript_id', 'gene_id']],
            on='transcript_id',
            how='left'
        )
        final_output.rename(columns={'Synt_id_x': 'Synt_id'}, inplace=True)
        # print duplicated rows
        print(final_output[final_output.duplicated()])
        # drop duplicated rows
        # final_output = final_output.drop_duplicates()
        # print rows where transcript_id is duplicated
        print(final_output[final_output.index.duplicated()])
        # set gene_id as index
        final_output.set_index('gene_id', inplace=True)
        # Save to TSV file
        final_output.to_csv(f'{args.output}_categories.tsv', sep='\t', index=True)

    print("Analysis complete!")


if __name__ == '__main__':
    main()



# python /DKED/scratch/nadjafn/potato-allelic-orthogroups/scripts/parse_pangenes.py --pangenes De_v1.functional_annotations_nodigits_genespace.tsv                                               --gff De_v1.functional_annotations_nodigits.longest_isoforms.gff                                                          --output De_v1.functional_annotations_nodigits_genespace -s 1hap1_1hap2_1hap3_1hap4_synteny
