import pandas as pd
import argparse
import re

def extract_number_before_G(s):
    # Regular expression pattern to find digits followed by 'G'
    pattern = r'(\d+)G'
    match = re.search(pattern, s)
    if match:
        return match.group(1)  # Return the matched digits
    else:
        return None  # No match found

def read_file(file_path):
    with open(file_path, 'r') as file:
        return file.read()


def parse_file_to_dataframe(file_content):
    # Initialize an empty list to store the data
    data = []

    # Current cluster number
    current_cluster = None

    # Iterate through each line in the file content
    for line in file_content.split('\n'):
        # Check if the line indicates a new cluster
        if line.startswith('>Cluster'):
            # Extract the cluster number
            current_cluster = int(line.split(' ')[1])
        else:
            # Extract the gene ID (if the line is not empty)
            if line.strip():
                gene_id = line.split('>')[1].split()[0]
                gene_id = gene_id.replace('...', '')
                # Extracting numbers
                haplotype = extract_number_before_G(gene_id)
                # Append the cluster number and gene ID to the data list
                data.append({'cluster': current_cluster, 'gene_id': gene_id, 'haplotype': haplotype})

    # Create a DataFrame from the data
    df = pd.DataFrame(data)
    return df

def analyze_dataframe(df):
    # Calculate the number of genes in each cluster
    cluster_counts = df['cluster'].value_counts().rename_axis('cluster').reset_index(name='num_genes')
    print(cluster_counts)
    # calculate the number of clusters with one to more genes
    in_cluster_counts = cluster_counts['num_genes'].value_counts().rename_axis('genes_in_cluster').reset_index(name='num_clusters')
    print(in_cluster_counts)
    # Group by 'cluster' and calculate the required statistics
    cluster_stats = df.groupby('cluster').agg(
        num_genes=('gene_id', 'count'), 
        num_haplotypes=('haplotype', 'nunique')
    ).reset_index()
    print(cluster_stats)

    # Calculate the combinations of number of genes and number of haplotypes in each cluster
    combination_counts = cluster_stats.groupby(['num_genes', 'num_haplotypes']).size().reset_index(name='count')
    print(combination_counts)

def main():
    parser = argparse.ArgumentParser(description="Parse a file to create a DataFrame.")
    parser.add_argument("file_path", help="Path to the file to be parsed")
    args = parser.parse_args()

    file_content = read_file(args.file_path)
    df = parse_file_to_dataframe(file_content)
    print(df)
    analyze_dataframe(df)

if __name__ == "__main__":
    main()