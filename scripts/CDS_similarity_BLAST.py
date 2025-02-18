import pandas as pd
import argparse
import numpy as np
import matplotlib.pyplot as plt

def parse_blast_fmt6(blast_file):
    """
    Parse the blast file in format 6
    :param blast_file: the blast file in format 6
    :return: a pandas dataframe
    """
    blast_df = pd.read_csv(blast_file, sep='\t', header=None)
    blast_df.columns = ['query', 'subject', 'identity', 'length', 'mismatch', 'gap', 'q_start', 'q_end', 's_start',
                        's_end', 'evalue', 'bitscore']
    return blast_df

def filter_blast(blast_df):
    # remove rows where query and subject are the same
    blast_df = blast_df[blast_df['query'] != blast_df['subject']]

    

    return blast_df



if __name__ == '__main__':
    # argparser to parse the input file

    parser = argparse.ArgumentParser(description='Parse the pagenome file and gff file')

    parser.add_argument('-b', '--blast', help='The blast all-vs-all', required=True)
    parser.add_argument('-s', '--synthenic_genes', help='Tsv file with syntenic genes', required=True)
    #parser.add_argument('-o', '--output', help='output_prefix', required=True)

    args = parser.parse_args()

    # parse the pan-genome file
    blast = parse_blast_fmt6(args.blast)
    #print(blast.head())
    filtered_blast = filter_blast(blast)


    syntelogs = pd.read_csv(args.synthenic_genes, sep='\t', usecols=['Synt_id', 'haplotype_id', 'CDS_haplotype_with_longest_annotation'])

    syntelogs = syntelogs[syntelogs['CDS_haplotype_with_longest_annotation'] == 'equal_lengths']

    syntelogs['haplotype'] = syntelogs['haplotype_id'].str.extract(r'(H\d+)')

    # group syntelogs by Synt_id
    syntelogs_grouped = syntelogs.groupby('Synt_id')['haplotype_id'].apply(lambda x: ','.join(x)).reset_index()
    

    # merge the blast and syntelogs
    merged = pd.merge(filtered_blast, syntelogs, left_on='query', right_on='haplotype_id', how = 'inner')
    merged2 = pd.merge(merged, syntelogs, left_on='subject', right_on='haplotype_id', how = 'inner')

    # only keep where synt_id_x == synt_id_y
    merged2 = merged2[merged2['Synt_id_x'] == merged2['Synt_id_y']]


    # get the number of unique synt_idds
    print(len(merged2['Synt_id_x'].unique()))
    # only keep haplotype_x and haplotype_y Synt_id_x and identiy and match
    
    # make a new column with a combination of haplotype_x and haplotype_y ordered alphabetically
    merged2['haplotype_comb'] = merged2[['haplotype_x', 'haplotype_y']].apply(lambda x: '-'.join(sorted(x)), axis=1)
 
    merged2.to_csv('results/03_GENESPACE/De_v1.functional_annotations_nodigits_genespace_syntelogs_all_exploded_blast.tsv', sep='\t', index=False)
    merged2 = merged2[['haplotype_comb', 'Synt_id_x', 'identity', 'mismatch']]

    # sort by identity and mismatch
    merged2 = merged2.sort_values(by=['identity', 'mismatch'])
    # drop duplicates
    merged2 = merged2.drop_duplicates(['haplotype_comb', 'Synt_id_x'], keep='first')
    # count per synt id how many of the haplotype combinations are identical
    merged2 = merged2.groupby(['Synt_id_x', 'identity', 'mismatch']).size().reset_index(name='counts_of_identical_combinations')

    # if identy is 100% and mismatch is 0, then the haplotypes are identical

    merged2['mismatch_category'] = np.where((merged2['identity'] == 100.000) & (merged2['mismatch'] == 0), 'identical', 'SNPs')
# for each Synt_id_x where mismatch_category is identical is not present then add a row with 0 counts
    synt_ids = merged2['Synt_id_x'].unique()
    for synt_id in synt_ids:
        if 'identical' not in merged2[merged2['Synt_id_x'] == synt_id]['mismatch_category'].values:
            new_row = pd.DataFrame([{'Synt_id_x': synt_id, 'identity': 100.000, 'mismatch': 0, 'counts_of_identical_combinations': 0, 'mismatch_category': 'identical'}])
            merged2 = pd.concat([merged2, new_row], ignore_index=True)
    # sort on Synt_id
    merged2 = merged2.sort_values(by='Synt_id_x')
    
    # Create a pivot table for plotting
    merged2_ident = merged2[merged2['mismatch_category'] == 'identical']
    # group by counts and count the occurence
    merged2_ident = merged2_ident.groupby(['counts_of_identical_combinations','mismatch']).size().reset_index(name='counts_per_synt_id')


    # add a category name for the counts
    merged2_ident['category'] = 'no_cat'
    merged2_ident['category'] = np.where((merged2_ident['counts_of_identical_combinations'] == 0), 'all_alleles_different', merged2_ident['category'])
    merged2_ident['category'] = np.where((merged2_ident['counts_of_identical_combinations'] == 1), '2_alleles_identical', merged2_ident['category'])
    merged2_ident['category'] = np.where((merged2_ident['counts_of_identical_combinations'] == 2), '2x2_alleles_identical', merged2_ident['category'])
    merged2_ident['category'] = np.where((merged2_ident['counts_of_identical_combinations'] == 3), '3_alleles_identical', merged2_ident['category'])
    merged2_ident['category'] = np.where((merged2_ident['counts_of_identical_combinations'] == 6), '4_alleles_identical', merged2_ident['category'])
    print(merged2_ident)

    # make a barplot of the counts
    merged2_ident.plot(kind='bar', x='category', y='counts_per_synt_id')
    # Plot the histogram
    # merged2.pivot(index='', columns='mismatch_category', values='counts').plot(kind='bar', stacked=True)
    # plt.xlabel('Syntenic ID')
    # plt.ylabel('Counts')
    # plt.title('Number of Identical and SNP Haplotypes per Syntenic ID')
    # plt.legend(title='Mismatch Category')
    plt.tight_layout()

    # # Save the plot
    plt.savefig('results/03_GENESPACE/haplotype_histogram.png')
    plt.close()

    # print 
    print(merged2)
    merged2_SNPs = merged2[merged2['mismatch_category'] == 'SNPs']
    print(merged2_SNPs)
    # make a histogram of the SNPs
    merged2_SNPs['mismatch'] = merged2_SNPs['mismatch'].astype(int)
    counts_mismatches = merged2_SNPs['mismatch'].value_counts()
    print(counts_mismatches)
    counts_mismatches.columns = ['mismatch', 'counts']
    counts_mismatches.plot(kind='hist', x='mismatch', y='counts', bins=100)
    plt.xlabel('Number of SNPs')
    plt.ylabel('Counts')
    plt.title('Number of SNPs per Syntenic ID')

    plt.tight_layout()
    plt.savefig('results/03_GENESPACE/SNP_histogram.png')


    # save the merged file


# python /scratch/nadjafn/potato-allelic-orthogroups/scripts/CDS_similarity_BLAST.py  -b /scratch/nadjafn/potato-allelic-orthogroups/work/98/25fb94d5a4ffc0b2f4e7e107477501/De_v1.functional_annotations_nodigits.txt  -s results/03_GENESPACE/De_v1.functional_annotations_nodigits_genespace_syntelogs_all_exploded.tsv 