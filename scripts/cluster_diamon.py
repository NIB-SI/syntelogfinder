import argparse
import pandas as pd
import networkx as nx
import re
import matplotlib.pyplot as plt


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
    
def create_clusters(input_file, output_file):
    # Read the TSV file
    df = pd.read_csv(input_file, sep="\t", header=None, usecols=[0, 1])

    # Create a graph from the DataFrame
    G = nx.from_pandas_edgelist(df, source=0, target=1)

    # Find connected components (clusters)
    clusters = list(nx.connected_components(G))

    # Create a new DataFrame for output
    cluster_data = []
    for cluster_id, genes in enumerate(clusters):
        for gene in genes:
            haplotype = extract_number_before_G(gene)
            cluster_data.append({'cluster': cluster_id, 'gene_id': gene, 'haplotype': haplotype})

    cluster_df = pd.DataFrame(cluster_data)

    # Write to a new file
    cluster_df.to_csv(output_file, index=False, sep='\t')
    return cluster_df

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
    return combination_counts


def plot_bar(df):
    # Pivoting data for stacked bar plot
    pivot_df = df.pivot(index='num_genes', columns='num_haplotypes', values='count').fillna(0)

    # Normal scale plot with increased label sizes
    ax1 = pivot_df.plot(kind='bar', stacked=True, figsize=(10, 6))
    ax1.set_xlabel('Number of Genes in cluster', fontsize=12)
    ax1.set_ylabel('Count', fontsize=12)
    ax1.set_title('Normal Scale: Gene Count vs. Number of Haplotypes in Cluster', fontsize=14)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.savefig('nolog.png')
    plt.close()  # Close the plot to prevent it from displaying here

    # Log scale plot with increased label sizes
    ax2 = pivot_df.plot(kind='bar', stacked=True, figsize=(10, 6), logy=True)
    ax2.set_xlabel('Number of Genes in cluster', fontsize=12)
    ax2.set_ylabel('Count (log scale)', fontsize=12)
    ax2.set_title('Log Scale: Gene Count vs. Number of Haplotypes in Cluster', fontsize=14)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.savefig('log.png')


def main():
    parser = argparse.ArgumentParser(description="Cluster genes based on their mappings.")
    parser.add_argument("input_file", help="Path to the input TSV file")
    parser.add_argument("output_file", help="Path to the output file")
    args = parser.parse_args()

    df = create_clusters(args.input_file, args.output_file)
    counts = analyze_dataframe(df)
    plot_bar(counts)

if __name__ == "__main__":
    main()
