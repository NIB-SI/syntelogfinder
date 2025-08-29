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

# Global constants
GENE_PATTERN = r'(\S+)$'
NA_VALUES = ['NA', 'N/A', '', None]
CUSTOM_ORDER = [
    'not_annotated', 'less_1%_difference', 'more_1%_difference',
    'more_5%_difference', 'more_10%_difference', 'more_20%_difference'
]

def parse_pangenes(pangenes_file: str) -> pd.DataFrame:
    """Parse the tab-separated pangenes file."""
    df = pd.read_csv(pangenes_file, sep="\t", header=0, index_col=False, na_values=NA_VALUES)
    return df.replace(', ', ',', regex=True)

def identify_overlapping_rows(df: pd.DataFrame, na_values: Optional[List] = None) -> Dict[int, Set[int]]:
    """Identify groups of rows that have overlapping gene IDs."""
    if na_values is None:
        na_values = NA_VALUES

    # Map gene IDs to the rows they appear in
    gene_to_rows = defaultdict(set)

    # Build the gene-to-rows mapping
    for row_idx, row in df.iterrows():
        for col in df.columns:
            value = row[col]
            if pd.isna(value) or value in na_values:
                continue

            if isinstance(value, str):
                matches = re.findall(GENE_PATTERN, value)
                for gene_id in matches:
                    gene_to_rows[gene_id].add(row_idx)

    # Create a graph of connected rows
    row_connections = defaultdict(set)
    for gene, rows in gene_to_rows.items():
        if len(rows) > 1:
            for row in rows:
                row_connections[row].update(rows)

    # Find connected components using BFS
    components = {}
    visited = set()

    for start_row in row_connections:
        if start_row in visited:
            continue

        queue = deque([start_row])
        component = set()

        while queue:
            row = queue.popleft()
            if row in visited:
                continue

            visited.add(row)
            component.add(row)

            for connected_row in row_connections[row]:
                if connected_row not in visited:
                    queue.append(connected_row)

        for row in component:
            components[row] = component

    # Add singleton rows
    for row_idx in df.index:
        if row_idx not in components:
            components[row_idx] = {row_idx}

    return components

def merge_overlapping_rows(df: pd.DataFrame, na_values: Optional[List] = None,
                          print_operations: bool = True) -> pd.DataFrame:
    """Merge rows that have overlapping gene IDs into consolidated rows."""
    if na_values is None:
        na_values = NA_VALUES

    components = identify_overlapping_rows(df, na_values)

    # Unique components
    unique_components = {}
    for component in components.values():
        rep = min(component)
        if rep not in unique_components:
            unique_components[rep] = component

    if print_operations:
        num_groups = len(unique_components)
        num_affected_rows = sum(len(comp) for comp in unique_components.values())
        num_merged = num_affected_rows - num_groups
        print(f"Found {num_groups} groups of overlapping rows.")
        print(f"Merging {num_merged} rows into their respective groups.\n")

    result_data = []

    for rep, component in unique_components.items():
        if len(component) == 1:
            result_data.append(df.loc[rep].to_dict())
            continue

        if print_operations:
            print(f"Merging rows: {', '.join(map(str, sorted(component)))}")

        merged_row = {}
        for col in df.columns:
            gene_ids = set()
            for row_idx in component:
                value = df.at[row_idx, col]
                if pd.isna(value) or value in na_values:
                    continue

                if isinstance(value, str):
                    matches = re.findall(GENE_PATTERN, value)
                    gene_ids.update(matches)

            gene_ids.discard('')
            merged_row[col] = ','.join(sorted(gene_ids)) if gene_ids else pd.NA

        result_data.append(merged_row)

    result_df = pd.DataFrame(result_data, columns=df.columns)

    if print_operations:
        print(f"\nOriginal dataframe had {len(df)} rows.")
        print(f"New dataframe has {len(result_df)} rows.")
        print(f"Reduced by {len(df) - len(result_df)} rows.")

    return result_df

def remove_duplicates(df: pd.DataFrame) -> pd.DataFrame:
    """Remove duplicate gene IDs within cells."""
    def process_cell(cell):
        if not isinstance(cell, str):
            return cell

        gene_ids = cell.split(',')
        unique_genes = {}

        for gene in gene_ids:
            base_gene = re.sub(r'[\*\+]', '', gene)
            if base_gene not in unique_genes or '*' in gene or '+' in gene:
                unique_genes[base_gene] = gene

        return ','.join(unique_genes.values())

    return df.map(process_cell)

def parse_gff(gff_file: str) -> pd.DataFrame:
    """Parse GFF file and extract relevant information."""
    gff = pd.read_csv(gff_file, sep="\t", header=None, index_col=False, comment='#')
    gff.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']

    # Extract transcript_id
    gff['transcript_id'] = ''

    # For mRNA/transcript entries (use transcript_id attribute)
    if 'mRNA' in gff['type'].values:
        mask_mrna = gff['type'].str.contains('mRNA', na=False)
    else:
        mask_mrna = gff['type'].str.contains('transcript', na=False)

    gff.loc[mask_mrna, 'transcript_id'] = gff.loc[mask_mrna, 'attributes'].str.extract(r'ID=([^;]+)')[0]

    # For features (use Parent attribute)
    mask_features = gff['type'].str.contains('exon|CDS|five_prime_UTR|three_prime_UTR', na=False)
    gff.loc[mask_features, 'transcript_id'] = gff.loc[mask_features, 'attributes'].str.extract(r'Parent=([^;]+)')[0]

    # Extract gene_id
    gff['gene_id'] = np.where(
        mask_mrna,
        gff['attributes'].str.extract(r'gene(?:_id|ID)=([^;]+)')[0],
        ''
    )

    # Create mapping and propagate gene_id
    transcript_to_gene = gff[gff['gene_id'] != ''].set_index('transcript_id')['gene_id'].to_dict()
    gff['gene_id'] = gff['transcript_id'].map(transcript_to_gene).fillna('')

    return gff

def count_comma_separated_values(df: pd.DataFrame, column_name: str) -> pd.Series:
    """Count the number of comma-separated values in a column."""
    def count_values(cell):
        if pd.isna(cell):
            return 0
        return len(str(cell).split(','))

    return df[column_name].apply(count_values)

def get_attribute_lengths(gff_file: str, attribute: str) -> pd.DataFrame:
    """Extract feature lengths from a GFF file for a specific attribute type."""
    db_filename = 'gff.db'

    if not os.path.exists(db_filename):
        print(f"Creating GFF database: {db_filename}")
        db = gffutils.create_db(gff_file, dbfn=db_filename, keep_order=False,
                                merge_strategy="error", disable_infer_transcripts=True, verbose=True)
    else:
        print(f"Using existing GFF database: {db_filename}")

    db = gffutils.FeatureDB(db_filename, keep_order=True)
    transcript_info = {}

    # Get all transcript IDs
    try:
        all_transcripts = {feature.id for feature in db.features_of_type("mRNA")}
    except:
        all_transcripts = {feature.id for feature in db.features_of_type("transcript")}

    # Process features
    for feature in db.features_of_type(attribute):
        parent_id = feature.attributes.get('Parent', [None])[0]
        if parent_id is None:
            parent_id = feature.attributes.get('ID', [None])[0]

        length = feature.end - feature.start + 1

        if parent_id not in transcript_info:
            try:
                parent = db[parent_id]
                parent_start = parent.start
            except gffutils.FeatureNotFoundError:
                parent_start = None

            transcript_info[parent_id] = {
                f'{attribute}_ref_length': 0,
                f'{attribute}_parent_start': parent_start
            }

        transcript_info[parent_id][f'{attribute}_ref_length'] += length

    # Ensure all transcripts are included
    for transcript_id in all_transcripts:
        if transcript_id not in transcript_info:
            transcript_info[transcript_id] = {
                f'{attribute}_ref_length': 0,
                f'{attribute}_parent_start': None
            }

    df = pd.DataFrame.from_dict(transcript_info, orient='index')
    df.reset_index(inplace=True)
    df.rename(columns={'index': 'transcript_id'}, inplace=True)

    return df

def get_haplotype_columns(df: pd.DataFrame) -> List[str]:
    """Get a list of haplotype columns from the DataFrame."""
    return [col for col in df.columns if col.startswith('hap')]

def get_haplotype_labels(df: pd.DataFrame) -> List[str]:
    """Get haplotype labels based on column names."""
    columns = [col for col in df.columns if '_length_' in col and col.endswith('G')]
    labels = [col.split('_')[-1] for col in columns]
    return sorted(labels)

def add_synteny_category(pangenes: pd.DataFrame) -> pd.DataFrame:
    """Add synteny categories to the pangenes DataFrame."""
    haplotype_cols = get_haplotype_columns(pangenes)


    # Count genes in each haplotype
    for hap in haplotype_cols:
        pangenes[f'{hap}_count'] = count_comma_separated_values(pangenes, hap).astype(str) + hap

    # Add synteny classification
    pangenes['true_synteny'] = 's'

    # Check for special characters
    has_special_chars = False
    for hap in haplotype_cols:
        has_special_chars = has_special_chars | pangenes[hap].str.contains('\*') | pangenes[hap].str.contains('\+')

    pangenes.loc[has_special_chars, 'true_synteny'] = 'no_s'

    # Add synteny ID
    pangenes['Synt_id'] = 'Synt_id_' + pangenes.index.astype(str)

    # Create synteny category
    synteny_category_parts = [pangenes[f'{hap}_count'] for hap in haplotype_cols]
    synteny_category_parts = pd.DataFrame(synteny_category_parts).T
    pangenes['synteny_category'] = synteny_category_parts.astype(str).agg('_'.join, axis=1) + '_' + pangenes['true_synteny']

    # Create clean copy without special characters
    pangenes_gene = pangenes.map(lambda x: x.replace("+", "").replace("*", "") if isinstance(x, str) else x)

    # Combine all haplotype columns
    syntenic_genes_parts = [pangenes_gene[hap] for hap in haplotype_cols]
    syntenic_genes_parts = pd.DataFrame(syntenic_genes_parts).T
    pangenes_gene['syntenic_genes'] = syntenic_genes_parts.astype(str).agg(
        lambda x: ','.join(v for v in x if pd.notna(v) and v not in NA_VALUES), axis=1
    )


    # Pivot the table
    pangenes_pivot = pd.melt(
        pangenes_gene,
        id_vars=['Synt_id', 'synteny_category', 'syntenic_genes'],
        value_vars=haplotype_cols,
        var_name='haplotype',
        value_name='transcript_id'
    )

    # Clean and expand
    pangenes_pivot = pangenes_pivot.dropna(subset='transcript_id')
    pangenes_pivot['transcript_id'] = pangenes_pivot['transcript_id'].str.split(',')
    pangenes_pivot = pangenes_pivot.explode('transcript_id')
    pangenes_pivot = pangenes_pivot.dropna(subset=['transcript_id'])
    pangenes_pivot = pangenes_pivot[pangenes_pivot['transcript_id'].str.strip() != '']
    pangenes_pivot = pangenes_pivot.drop_duplicates(subset=['transcript_id'])

    return pangenes_pivot

def merge_pangenes_gff(pangenes_pivot: pd.DataFrame, gff: pd.DataFrame) -> pd.DataFrame:
    """Merge pangenes and GFF DataFrames."""
    # try frist mRNA, then transcript
    if not gff['type'].isin(['mRNA', 'transcript']).any():
        raise ValueError("GFF file does not contain 'mRNA' or 'transcript' entries.")
    if 'mRNA' in gff['type'].values:
        gff_mrna = gff[gff['type'] == 'mRNA'].copy()
    else:
        gff_mrna = gff[gff['type'] == ['transcript']].copy()

    # Ensure consistent data types
    gff_mrna['transcript_id'] = gff_mrna['transcript_id'].astype(str)
    pangenes_pivot['transcript_id'] = pangenes_pivot['transcript_id'].astype(str)

    gff_pangenes = pd.merge(gff_mrna, pangenes_pivot, on='transcript_id', how='left')
    gff_pangenes['synteny_category'] = gff_pangenes['synteny_category'].astype(str)

    return gff_pangenes

def percent_func(pct: float, allvalues: np.ndarray) -> str:
    """Format percentage labels for pie charts."""
    absolute = int(pct / 100. * np.sum(allvalues))
    return "{:.1f}%\n({:d})".format(pct, absolute)

def make_pie_chart(gff_pangenes: pd.DataFrame, syntelogs_category: str, output_prefix: str) -> None:
    """Create a pie chart showing distribution of synteny categories."""
    mrna_data = gff_pangenes[gff_pangenes['type'].isin(['mRNA', 'transcript'])].copy()
    mrna_data = mrna_data.sort_values('synteny_category')

    synt_counts = mrna_data['synteny_category'].value_counts()


    # Get top 7 categories and group the rest as "other"
    synt_counts_top = synt_counts[:7].copy()
    if len(synt_counts) > 7:
        synt_counts_top['other'] = synt_counts[7:].sum()
    synt_counts_top.sort_index(inplace=True)

    # Set up colors
    colors = ['#D3D3D3'] * len(synt_counts_top)
    if syntelogs_category in synt_counts_top.index:
        highlight_idx = list(synt_counts_top.index).index(syntelogs_category)
        colors[highlight_idx] = '#FF0000'

    # Create pie chart
    explode = tuple([0.1 * (7-i) for i in range(len(synt_counts_top))])
    plt.figure(figsize=(7, 7))
    synt_counts_top.plot.pie(
        startangle=90,
        explode=explode,
        autopct=lambda pct: percent_func(pct, synt_counts),
        colors=colors
    )

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_pie_chart.svg', bbox_inches='tight')
    plt.close()

def check_length_values(row: pd.Series, percent: float) -> bool:
    """Check if the difference between min and max values exceeds the specified percentage."""
    min_value = row.min()
    max_value = row.max()
    cutoff = min_value * percent * 0.01
    return max_value >= min_value + cutoff

def add_length_analysis(df: pd.DataFrame, attribute: str) -> pd.DataFrame:
    """Add length categories, differences, and longest transcript analysis."""
    # Select length columns
    length_cols = [col for col in df.columns if f'{attribute}_length_' in col and col.endswith('G')]
    df[length_cols] = df[length_cols].fillna(0)
    df_subset = df[length_cols].astype(int)

    # Initialize category column
    df[f'{attribute}_length_category'] = 'unclassified'

    # Apply categorization rules
    for threshold, category in [(1, 'less_1%_difference'), (5, 'more_1%_difference'),
                               (10, 'more_5%_difference'), (20, 'more_10%_difference')]:
        if category == 'less_1%_difference':
            mask = ~df_subset.apply(lambda x: check_length_values(x, 1), axis=1)
        elif category == 'more_1%_difference':
            mask = (df_subset.apply(lambda x: check_length_values(x, 1), axis=1) &
                   ~df_subset.apply(lambda x: check_length_values(x, 5), axis=1))
        elif category == 'more_5%_difference':
            mask = (df_subset.apply(lambda x: check_length_values(x, 5), axis=1) &
                   ~df_subset.apply(lambda x: check_length_values(x, 10), axis=1))
        elif category == 'more_10%_difference':
            mask = (df_subset.apply(lambda x: check_length_values(x, 10), axis=1) &
                   ~df_subset.apply(lambda x: check_length_values(x, 20), axis=1))

        df.loc[mask, f'{attribute}_length_category'] = category

    # Handle special cases
    df.loc[df_subset.apply(lambda x: check_length_values(x, 20), axis=1),
           f'{attribute}_length_category'] = 'more_20%_difference'
    df.loc[(df_subset == 0).all(axis=1), f'{attribute}_length_category'] = 'not_annotated'

    # Add percentage differences
    df[f'{attribute}_max_difference'] = df_subset.max(axis=1) - df_subset.min(axis=1)
    df[f'{attribute}_percent_difference'] = (df[f'{attribute}_max_difference'] /
                                           df_subset.min(axis=1) * 100)

    # Add longest transcript information
    df[f'{attribute}_haplotype_with_longest_annotation'] = df[length_cols].idxmax(axis=1)
    mask = (df[length_cols].nunique(axis=1) == 1)
    df.loc[mask, f'{attribute}_haplotype_with_longest_annotation'] = 'equal_lengths'
    df[f'{attribute}_haplotype_with_longest_annotation'] = (
        df[f'{attribute}_haplotype_with_longest_annotation'].str.replace(f'{attribute}_length_', '')
    )

    return df

def make_barplot(df: pd.DataFrame, attribute: str, ax: plt.Axes, haplotype_labels: List[str]) -> None:
    """Create a barplot for length categories on the given axes."""
    # Convert to categorical for consistent ordering
    df[f'{attribute}_length_category'] = pd.Categorical(
        df[f'{attribute}_length_category'],
        categories=CUSTOM_ORDER,
        ordered=True
    )

    haplotype_categories = haplotype_labels + ['equal_lengths']
    df[f'{attribute}_haplotype_with_longest_annotation'] = pd.Categorical(
        df[f'{attribute}_haplotype_with_longest_annotation'],
        categories=haplotype_categories,
        ordered=True
    )

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
    ax.tick_params(axis='x', rotation=90)

def generate_combined_plots(df_list: List[pd.DataFrame], attributes: List[str], output_prefix: str) -> None:
    """Generate combined barplots for multiple attributes."""
    fig, axes = plt.subplots(1, len(attributes), figsize=(13, 5), sharey=True)

    if len(attributes) == 1:
        axes = [axes]

    for i, (df, attribute) in enumerate(zip(df_list, attributes)):
        make_barplot(df, attribute, axes[i], get_haplotype_labels(df))
        if i > 0:
            axes[i].legend_.remove()

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_combined_barplots.svg')
    plt.close(fig)

def filter_all_syntelogs(pangenes: pd.DataFrame, target_category: str = None) -> pd.DataFrame:
    """Filter pangenes to include only specified synteny category."""
    if target_category is None:
        haplotype_cols = [col for col in pangenes.columns if col.startswith('hap')]
        filter_parts = [f"1{col}" for col in haplotype_cols]
        target_category = '_'.join(filter_parts) + '_s'

    return pangenes[pangenes['synteny_category'] == target_category].copy()

def add_length_to_syntelogs(pangenes: pd.DataFrame, ref_lengths: pd.DataFrame, attribute: str) -> pd.DataFrame:
    """Add length information to syntelog groups."""
    df_synt_lengths = pd.merge(pangenes, ref_lengths, on='transcript_id', how='inner')

    transcript_ids = (df_synt_lengths.groupby('Synt_id')['transcript_id']
                     .apply(list).reset_index())

    haplotypes = df_synt_lengths['haplotype'].unique()

    # Pivot data
    df_synt_pivot = df_synt_lengths.pivot(
        index='Synt_id',
        columns=['haplotype'],
        values=[f'{attribute}_ref_length', f'{attribute}_parent_start']
    )

    df_synt_pivot.columns = ['_'.join(col) for col in df_synt_pivot.columns]
    df_synt_pivot = pd.merge(df_synt_pivot, transcript_ids, on='Synt_id', how='inner')

    # Rename columns
    renamed_columns = ['Synt_id']
    for hap in sorted(haplotypes):
        hap_num = hap.replace('hap', '')
        renamed_columns.append(f'{attribute}_length_{hap_num}G')
    for hap in sorted(haplotypes):
        hap_num = hap.replace('hap', '')
        renamed_columns.append(f'{attribute}_parent_start_{hap_num}G')
    renamed_columns.append('transcript_id')

    df_synt_pivot.columns = renamed_columns

    # Add length analysis
    return add_length_analysis(df_synt_pivot, attribute)

def main():
    """Main function to execute the analysis workflow."""
    parser = argparse.ArgumentParser(description='Parse the pan-genome file and GFF file')
    parser.add_argument('-p', '--pangenes', help='The pan-genome file', required=True)
    parser.add_argument('-g', '--gff', help='The GFF file', required=True)
    parser.add_argument('-o', '--output', help='Output prefix', required=True)
    parser.add_argument('-s', '--syntelogs_category', help='Target syntelog category to filter for', default=None)

    args = parser.parse_args()

    # Parse input files
    print("Parsing pan-genome file...")
    pangenes = parse_pangenes(args.pangenes)

    print("Merging overlapping rows...")
    pangenes = merge_overlapping_rows(pangenes, print_operations=False)
    pangenes = remove_duplicates(pangenes)
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

    # Merge pangenes and GFF data
    print("Merging pangenes and GFF data...")
    gff_pangenes = merge_pangenes_gff(pangenes_pivot, gff)

    # Create pie charts
    print("Creating pie chart...")
    make_pie_chart(gff_pangenes, args.syntelogs_category, args.output)

    # Filter for syntelogs
    print("Filtering syntelogs...")
    syntelogs = filter_all_syntelogs(pangenes_pivot, args.syntelogs_category)

    # Process each attribute type
    print("Processing feature lengths...")
    attributes = ['exon', 'CDS']
    syntelogs_lengths_list = []

    for attribute in attributes:
        print(f"Processing {attribute} features...")
        ref_lengths = get_attribute_lengths(args.gff, attribute)
        syntelogs_lengths = add_length_to_syntelogs(syntelogs, ref_lengths, attribute)
        syntelogs_lengths_list.append(syntelogs_lengths)

    # Generate combined plots
    print("Generating combined plots...")
    generate_combined_plots(syntelogs_lengths_list, attributes, args.output)

    # Prepare final output
    print("Merging and saving final output...")
    if len(syntelogs_lengths_list) > 1:
        syntelogs_lengths_cds = syntelogs_lengths_list[1].explode('transcript_id')

        pangenes_pivot = pd.merge(
            pangenes_pivot,
            syntelogs_lengths_cds,
            on='transcript_id',
            how='left'
        )
        pangenes_pivot = pangenes_pivot.drop_duplicates('transcript_id')
        pangenes_pivot.set_index('transcript_id', inplace=True)

        final_output = pangenes_pivot[[
            'Synt_id_x', 'synteny_category', 'syntenic_genes', 'haplotype',
            'CDS_length_category', 'CDS_haplotype_with_longest_annotation',
            'exon_length_category', 'exon_haplotype_with_longest_annotation'
        ]]

        final_output = pd.merge(
            final_output,
            gff_pangenes[['transcript_id', 'gene_id']],
            on='transcript_id',
            how='left'
        )

        final_output.rename(columns={'Synt_id_x': 'Synt_id'}, inplace=True)
        final_output.set_index('gene_id', inplace=True)
        final_output.to_csv(f'{args.output}_categories.tsv', sep='\t', index=True)

    print("Analysis complete!")

if __name__ == '__main__':
    main()
