import pandas as pd
import argparse
import matplotlib.pyplot as plt


def clean_cell(cell):
    if pd.isna(cell):
        return cell  # Return as-is if the entire cell is NaN
    # Split the cell, filter out 'NaN', and join back into a string
    return ', '.join(filter(lambda x: x != 'NaN', cell.split(', ')))

def split_columns(df):
    # Columns that need to be split
    columns_to_split = ['Aa', 'Ab', 'Ac', 'Ad']

    for col in columns_to_split:
        # Splitting the column on ', ' and expanding the result into separate columns
        split_cols = df[col].str.split(', ', expand=True)

        # Naming the new columns
        split_cols.columns = [f"{col}.{i}" if i > 0 else col for i in range(split_cols.shape[1])]
        # Dropping the original column from the DataFrame
        df = df.drop(col, axis=1)

        # Joining the new columns to the DataFrame
        df = df.join(split_cols)

    return df 
    

# Function to concatenate columns with handling NaN values
def concatenate_columns(row, columns):
    values = []
    for col in columns:
        if not pd.isna(row[col]):
            values.append(str(row[col]))  # Convert to string
        else:
            values.append('NaN')  # Or any other string value you prefer for missing data


    return ', '.join(values)

# Function to count unique gene IDs per row
def count_unique_gene_ids(row):
    # Concatenating and splitting the gene IDs
    gene_ids = sum((str(row[col]).split(', ') for col in ['Aa', 'Ab', 'Ac', 'Ad']), [])
    # Removing empty strings and counting unique values
    return len(set(gene_id for gene_id in gene_ids if gene_id))

# Function to count the number of columns with gene IDs per row
def count_columns_with_gene_ids(row):
    return sum(pd.notnull(row[col]) and row[col] != '' for col in ['Aa', 'Ab', 'Ac', 'Ad'])


def filter_chromosomes(df):
    # Function to extract chromosome number
    def extract_chromosome(s):
        try:
            return s.split('_')[1].split('.')[1]
        except:
            return None

    # Apply the function to extract chromosome numbers
    df['Chromosome_Aa'] = df['Aa'].apply(extract_chromosome)
    df['Chromosome_Ab'] = df['Ab'].apply(extract_chromosome)
    df['Chromosome_Ac'] = df['Ac'].apply(extract_chromosome)
    df['Chromosome_Ad'] = df['Ad'].apply(extract_chromosome)
    
    # Add a column for the unique number of non-None chromosomes
    df['Unique_Chromosomes'] = df[['Chromosome_Aa', 'Chromosome_Ab', 'Chromosome_Ac', 'Chromosome_Ad']].apply(lambda row: len(set(row) - {None}), axis=1)


    # Helper function to compare two chromosomes, considering None as a wildcard
    def compare_chromosomes(chrom1, chrom2):
        return chrom1 == chrom2 or chrom1 is None or chrom2 is None

    # # Filter rows based on the new comparison logic
    # filtered_df = df[
    #     df.apply(lambda x: compare_chromosomes(x['Chromosome_Aa'], x['Chromosome_Ab']) and
    #                       compare_chromosomes(x['Chromosome_Ac'], x['Chromosome_Ad']) and
    #                       compare_chromosomes(x['Chromosome_Aa'], x['Chromosome_Ac']), axis=1)
    # ]
    # # Drop the temporary columns
    # filtered_df = filtered_df.drop(['Chromosome_Aa','Chromosome_Ab', 'Chromosome_Ac', 'Chromosome_Ad'], axis=1)
    df = df.drop(['Chromosome_Aa','Chromosome_Ab', 'Chromosome_Ac', 'Chromosome_Ad'], axis=1)

    return df

def merge_ortho_with_dup(df_orth, df_cluster):
    
    df_orth = split_columns(df_orth)

    # Melting the main dataframe to have all gene IDs in one column
    df_orth_melted = df_orth.melt(value_name="gene_id", var_name="column", ignore_index=False)

    # Performing a left join with the cluster dataframe
    merged_df = df_orth_melted.merge(df_cluster, on="gene_id", how="left")

    # Pivoting the result back to the original format
    # Adding cluster and haplotype columns to each original column
    pivot_columns = ['Aa', 'Aa.1', 'Aa.2', 'Aa.3', 'Aa.4', 'Ab', 'Ab.1', 'Ab.2', 'Ac', 'Ac.1', 'Ac.2', 'Ad']
    result_df = pd.DataFrame(index=df_orth.index)
    
    for col in pivot_columns:
        temp_df = merged_df[merged_df['column'] == col]
        result_df[col] = temp_df['gene_id'].values
        result_df[col + "_cluster"] = temp_df['cluster'].values
   
    # Apply the merging process for cluster columns
    for base_col in ['Aa', 'Ab', 'Ac', 'Ad']:
        # Identify all cluster columns that are variations of the base column (like 'Aa_cluster', 'Aa.1_cluster', ...)
        gene_id_columns = [col for col in result_df.columns if col.startswith(base_col) and 'cluster' not in col and 'haplotype' not in col]
        cluster_columns = [col for col in result_df.columns if col.startswith(base_col) and 'cluster' in col]

        # Create a new cluster column by concatenating all related cluster columns
        if cluster_columns:
            result_df[base_col + "_cluster"] = result_df.apply(lambda row: concatenate_columns(row, cluster_columns), axis=1)
            # Drop the individual cluster columns
            result_df = result_df.drop(columns=cluster_columns[1:])

          # Create a new cluster column by concatenating all related cluster columns
        if gene_id_columns:
            result_df[base_col] = result_df.apply(lambda row: concatenate_columns(row, gene_id_columns), axis=1)
            # Drop the individual cluster columns
            result_df = result_df.drop(columns=gene_id_columns[1:])

    # Apply the cleaning function to each cell in the DataFrame
    for col in result_df.columns:
        result_df[col] = result_df[col].apply(clean_cell)


    
    return result_df
   



def count_dups(result_df):
    # Select the cluster columns
    cluster_columns = ['Aa_cluster', 'Ab_cluster', 'Ac_cluster', 'Ad_cluster']

    # Count unique non-NaN values in each row for the cluster columns
    #print(result_df['Aa_cluster'][1].dropna().unique())
    # result_df['Unique_Sequences'] = result_df[cluster_columns].apply(lambda x: len(x.dropna().unique()), axis=1)

    result_df['Unique_Sequences'] = result_df[cluster_columns].apply(
    lambda row: len(set(
        item.strip() 
        for col in row.index 
        for item in str(row[col]).split(',') 
        if item.strip() != '' and not pd.isna(item)
    )), axis=1)



    result_df = filter_chromosomes(result_df)

    # Applying the functions to the dataframe
    result_df['Total_Gene_IDs'] = result_df.apply(count_unique_gene_ids, axis=1)
    result_df['Columns_with_Gene_IDs'] = result_df.apply(count_columns_with_gene_ids, axis=1)

    # Count the frequency of each unique cluster count
    cluster_count_freq = result_df['Unique_Sequences'].value_counts()

    # Create a stacked bar plot
    cluster_count_freq.plot(kind='bar', stacked=True)

    # Adding titles and labels
    plt.title('Distribution of Unique Sequences')
    plt.xlabel('Unique Cluster Count')
    plt.ylabel('Frequency')

    # Save the plot as a PNG file
    plt.savefig('unique_Sequence_counts.png')

    # Display the first few rows of the dataframe to verify
    result_df.to_csv('N0_short_dup.tsv', index=False, sep = "\t")


def get_perfect_homologs(data):
    filtered_data = data[~data.astype(str).apply(lambda row: row.str.contains(',')).any(axis=1)]
    print('# of lines without any gene dublications:', len(filtered_data))
    # Filter rows without any value in any column
    filtered_data = filtered_data.dropna(how='any')

    filtered = filter_chromosomes(filtered_data)
    return filtered
    
 
def read_duplication_file(filename):
    dup_file = pd.read_csv(filename, sep = "\t")
    return dup_file

def make_allelic_groups(filename):
    data = pd.read_csv(filename, sep = "\t")
    # Get orthogroups that have only one gene in them
    condition_1_gene = data.apply(lambda row: row.str.startswith('Soltu.Atl')).sum(axis=1) == 1
    one_gene_orthogroups = data.loc[condition_1_gene, data.columns[data.columns.str.startswith('Soltu.Atl') | (data.columns == 'OG')]]


    perfect_homolog = get_perfect_homologs(data)

    return perfect_homolog

def filter_csv(input_file, output_file):
    # Read the CSV file
    data = pd.read_csv(input_file, sep = "\t")
    # Create a condition to filter rows
    condition = data.apply(lambda row: row.str.startswith('Soltu.Atl')).sum(axis=1) == 4
    # Extract rows and select specific columns
    result_df = data.loc[condition, data.columns[data.columns.str.startswith('Soltu.Atl') | (data.columns == 'Orthogroup')]]

   
    
    # Filter rows that don't contain ',' or 'NA' in any column
    filtered_data = data[~data.astype(str).apply(lambda row: row.str.contains(',')).any(axis=1)]
    # Filter rows without any value in any column
    filtered_data = filtered_data.dropna(how='any')
    
    # Write the filtered data to a new CSV file
    filtered_data.to_csv(output_file, index=False,sep = "\t")
    

# Set up command-line arguments
parser = argparse.ArgumentParser(description='Filter a CSV file to remove rows containing "," or "NA".')
parser.add_argument('orthogroup_file', help='Input CSV file')
parser.add_argument('dublication_file', help='Output filtered CSV file')

args = parser.parse_args()
input_file = args.orthogroup_file
dup_file = args.dublication_file

# Filter the CSV file
# filter_csv(input_file, output_file)

#filtered_orthogroup = make_allelic_groups(input_file
df_dup = read_duplication_file(dup_file)
orthogroup_file = pd.read_csv(input_file, sep = "\t")
df = merge_ortho_with_dup(orthogroup_file, df_dup)
count_dups(df)
