# load packages
import pandas as pd
import re
import gffutils
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import argparse
import os


# build a function to parse the pan-genome file
def parse_pangenes(pangenes_file):
    return pd.read_csv(pangenes_file, sep="\t", header=0, index_col=False)

def parse_gff(gff_file):
    gff = pd.read_csv(gff_file, sep="\t", header=None, index_col=False, comment='#')
    gff.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    gff['ID'] = gff['attributes'].str.extract(r'ID=(.*?);')    

    return gff

def func(pct, allvalues):
        absolute = int(pct / 100.*np.sum(allvalues))
        return "{:.1f}%\n({:d})".format(pct, absolute)

def count_comma_separated_values(df, column_name):
    def count_values(cell):
        if pd.isna(cell):
            return 0
        return len(str(cell).split(','))

    return df[column_name].apply(count_values)


def get_transcript_lengths(gff_file, attribute):
    # Check if the gff.db file exists
    if not os.path.exists('gff.db'):
        # Create a database from the GFF file
        db = gffutils.create_db(gff_file, dbfn='gff.db', keep_order=False,
                                merge_strategy='merge', sort_attribute_values=False)
    
    # Load the existing gff.db file
    db = gffutils.FeatureDB('gff.db', keep_order=True)

    # Initialize a dictionary to store transcript lengths and start coordinates
    transcript_info = {}


    # Iterate over all features of type 'exon' in the GFF file
    for exon in db.features_of_type(attribute):
        parent_id = exon.attributes['Parent'][0]
        length = exon.end - exon.start + 1
        
        if parent_id not in transcript_info:
            # Get the parent transcript feature
            parent = db[parent_id]
            transcript_info[parent_id] = {f'{attribute}_ref_length': 0, f'{attribute}_parent_start': parent.start}
        
        transcript_info[parent_id][f'{attribute}_ref_length'] += length
    
    # Convert dictionary to DataFrame
    df = pd.DataFrame.from_dict(transcript_info, orient='index')
    df['haplotype_id'] = df.index
    return df

def make_barplot(df, attribute, output_prefix):
    plt.figure(figsize=(5, 5))
    # Sort the DataFrame by the 'length_category' column
    custom_order = ['less_1%_difference','more_1%_difference', 'more_5%_difference', 'more_10%_difference', 'more_20%_difference']  # replace with your actual categories
    df[f'{attribute}_length_category'] = pd.Categorical(df[f'{attribute}_length_category'], categories=custom_order, ordered=True)
    # chagne the haplotype_with_longest_annotation to a categorical variable
    df[f'{attribute}_haplotype_with_longest_annotation'] = pd.Categorical(df[f'{attribute}_haplotype_with_longest_annotation'], categories=['2G', '4G', '1G', '3G', 'equal_lengths'], ordered=True)
    #sns.countplot(x='length_category', hue='haplotype_with_longest_annotation', data=df)
    sns.histplot(data=df, x=f'{attribute}_length_category', hue=f'{attribute}_haplotype_with_longest_annotation', multiple='stack')
    plt.xlabel(f'{attribute} Length Category')
    plt.ylabel('Count')
    # turn the x-axis labels
    plt.xticks(rotation=90)
    plt.ylim(0, 15000)
    plt.tight_layout()
    # save the barplot
    plt.savefig(f'{output_prefix}_{attribute}_length_category_barplot.png')

    plt.show()

def add_synteny_category(pangenes):
    for hap in ['hap1', 'hap2', 'hap3', 'hap4']:
        pangenes[f'{hap}_count'] = count_comma_separated_values(pangenes, hap).astype(str) + hap
    # add a category if all the haplotypes dont't contain a * or + in the gene nam
    pangenes['true_synteny'] = 'synteny'
    pangenes.loc[(pangenes['hap1'].str.contains('\*') | pangenes['hap2'].str.contains('\*') | pangenes['hap3'].str.contains('\*') | pangenes['hap4'].str.contains('\*')), 'true_synteny'] = 'no_synteny'
    pangenes.loc[(pangenes['hap1'].str.contains('\+') | pangenes['hap2'].str.contains('\+') | pangenes['hap3'].str.contains('\+') | pangenes['hap4'].str.contains('\+')), 'true_synteny'] = 'no_synteny'

        # add a synt_id colum to the dataframe
    pangenes.reset_index(inplace=False)
    pangenes['Synt_id'] = pangenes.index
    # merge all the haplotype counts into one column
    pangenes['synteny_category'] = pangenes['hap1_count'] + '_' + pangenes['hap2_count'] + '_' + pangenes['hap3_count'] + '_' + pangenes['hap4_count'] + '_' + pangenes['true_synteny']
    pangenes_gene = pangenes.map(lambda x: x.replace("+", "").replace("*", "") if isinstance(x, str) else x)
    # print 
    pangenes_gene['syntenic_genes'] = pangenes_gene['hap1'] + ',' + pangenes_gene['hap2'] + ',' + pangenes_gene['hap3'] + ',' + pangenes_gene['hap4']
    # pivot the table
    pangenes_pivot = pd.melt(pangenes_gene, id_vars= ['Synt_id', 'synteny_category', 'syntenic_genes'], value_vars=['hap1', 'hap2', 'hap3', 'hap4'], var_name='haplotype', value_name='ID')
    # drop the rows where gene is NaN
    pangenes_pivot = pangenes_pivot.dropna(subset=['ID'])
    # convert gene column to a list
    pangenes_pivot['ID'] = pangenes_pivot['ID'].str.split(',')
    pangenes_pivot= pangenes_pivot.explode('ID')

    return pangenes_pivot

def merge_pangenes_gff(pangenes_pivot, gff):
    # select only mRNA from gff
    gff = gff[gff['type'] == 'mRNA']
    # Ensure the ID columns are of the same type (string)
    gff['ID'] = gff['ID'].astype(str)
    pangenes_pivot['ID'] = pangenes_pivot['ID'].astype(str)
    
    gff_pangenes = pd.merge(gff, pangenes_pivot, on='ID', how='left')
    # convert to string 
    gff_pangenes['synteny_category'] = gff_pangenes['synteny_category'].astype(str)
    return gff_pangenes

def make_pie_chart(gff_pangenes, output_prefix):
    # only select the rows that have mRNA in the type column
    gff_pangenes = gff_pangenes[gff_pangenes['type'] == 'mRNA']
    # sort the dataframe by the synteny category
    gff_pangenes = gff_pangenes.sort_values('synteny_category')
    # count the number of genes in each synteny category
    synt_counts = gff_pangenes['synteny_category'].value_counts()
    # if a gene is not in the top ten categories, replace it by "other"
    synt_counts_top = synt_counts[:5]
    synt_counts_top['other'] = synt_counts[5:].sum()
    # sort the dataframe by the synteny category alphabetically
    synt_counts_top.sort_index(inplace=True)

    explode = (0.6,0.4,0.3,0.1,0.2,0)
    
    # make a pie charr to show the distribution of synteny categories
    synt_counts_top.plot.pie(startangle=90, explode=explode, autopct= lambda pct: func(pct, synt_counts))
    # save the pie chart
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_pie_chart.png' )
    
def check_length_values(row, percent):
    min_value = row.min()
    max_value = row.max()

    # Calculate the 1% range
    cutoff_percent = min_value * percent * 0.01

    # Check if all values are within 1% of each other
    return max_value >= min_value + cutoff_percent


def add_length_category(df, attribute):
    # Select only the relevant columns
    df_subset = df[[f'{attribute}_length_1G', f'{attribute}_length_2G', f'{attribute}_length_3G', f'{attribute}_length_4G']]

    df_subset = df_subset.astype(int)

    # print max and min values of each row
    df[f'{attribute}_length_category'] = 'unclassified'

    df.loc[df_subset.apply(lambda x: check_length_values(x, 0), axis=1), f'{attribute}_length_category'] = 'less_1%_difference'
    df.loc[df_subset.apply(lambda x: check_length_values(x, 1), axis=1), f'{attribute}_length_category'] = 'more_1%_difference'

    df.loc[df_subset.apply(lambda x: check_length_values(x, 5), axis=1), f'{attribute}_length_category'] = 'more_5%_difference'

    df.loc[df_subset.apply(lambda x: check_length_values(x, 10), axis=1), f'{attribute}_length_category'] = 'more_10%_difference'


    df.loc[df_subset.apply(lambda x: check_length_values(x, 20), axis=1), f'{attribute}_length_category'] = 'more_20%_difference'
    return df

def add_length_differences_percent(df, attribute):
    df_subset = df[[f'{attribute}_length_1G', f'{attribute}_length_2G', f'{attribute}_length_3G', f'{attribute}_length_4G']]
    df_subset = df_subset.astype(int)

    # Calculate the differences between the columns
    # max difference
    df[f'{attribute}_max_difference'] = df_subset.max(axis=1) - df_subset.min(axis=1)

    # calculate the percentage difference
    df[f'{attribute}_percent_difference'] = df[f'{attribute}_max_difference'] / df_subset.min(axis=1) * 100

    return df

def add_longest_transcript(df, attribute):

    df[f'{attribute}_haplotype_with_longest_annotation'] = df[[f'{attribute}_length_1G', f'{attribute}_length_2G', f'{attribute}_length_3G', f'{attribute}_length_4G']].idxmax(axis=1)
    mask = (df[[f'{attribute}_length_1G', f'{attribute}_length_2G', f'{attribute}_length_3G', f'{attribute}_length_4G']].nunique(axis=1) == 1)
    df.loc[mask, f'{attribute}_haplotype_with_longest_annotation'] = 'equal_lengths'
    df[f'{attribute}_haplotype_with_longest_annotation'] = df[f'{attribute}_haplotype_with_longest_annotation'].str.replace(f'{attribute}_length_', '')
    #df.loc[df['length_category'] == 'less_1%_difference', 'haplotype_with_longest_annotation'] = 'lengths_within_5%'
    return df

def filter_4x_syntelogs(pangenes, filter = '1hap1_1hap2_1hap3_1hap4_synteny'):
    # select only the rows that have the category 1hap1_1hap2_1hap3_1hap4
    pangenes = pangenes[pangenes['synteny_category'] == filter]
    pangenes['haplotype_id'] = pangenes['ID']
    del pangenes['ID']
    return pangenes


def add_length_to_syntelogs(pangenes, ref_lengths, attribute):
    # Merge the two dataframes
    df_synt_lengths = pd.merge(pangenes, ref_lengths, on='haplotype_id', how='inner')

    # group the haplotypes ids by synteny id
    haplotype_ids = df_synt_lengths.groupby('Synt_id')['haplotype_id'].apply(list).reset_index()

    # sort by Synt_id and print duplicated ones
    haplotype_ids = haplotype_ids.sort_values(by='Synt_id')

    # print duplicated Synt_ids to check if there are any
    df_synt_pivot = df_synt_lengths.pivot(index='Synt_id', columns=['haplotype'], values=[f'{attribute}_ref_length',f'{attribute}_parent_start'])
    # flatten the columns
    df_synt_pivot.columns = ['_'.join(col) for col in df_synt_pivot.columns]
    df_synt_pivot = pd.merge(df_synt_pivot, haplotype_ids, on='Synt_id', how='inner')


    # # drop the synteny id column
    # df_synt_pivot.index = df_synt_pivot['Synt_id']
    # df_synt_pivot = df_synt_pivot.drop(columns=['Synt_id'])
    # # rename the columns
    # df_synt_pivot
    df_synt_pivot.columns = ['Synt_id',f'{attribute}_length_1G', f'{attribute}_length_2G', f'{attribute}_length_3G', f'{attribute}_length_4G',f'{attribute}_parent_start_1G', f'{attribute}_parent_start_2G', f'{attribute}_parent_start_3G', f'{attribute}_parent_start_4G', 'haplotype_id']
    # drop rows with missing values
    df_synt_pivot = df_synt_pivot.dropna()
    df_length = add_length_category(df_synt_pivot, attribute)
    df_length = add_length_differences_percent(df_length, attribute)

    df_length_cat = add_longest_transcript(df_length, attribute)

    return df_length_cat

def combine_attributes(row):
    attr_x = row['attributes_x'] if pd.notna(row['attributes_x']) else ''
    new_attr = row['new_attributes'] if pd.notna(row['new_attributes']) else ''
    return attr_x + new_attr






if __name__ == '__main__':
    # argparser to parse the input file

    parser = argparse.ArgumentParser(description='Parse the pagenome file and gff file')

    parser.add_argument('-p', '--pangenes', help='The pan-genome file', required=True)
    parser.add_argument('-g', '--gff', help='The gff file', required=True)
    parser.add_argument('-o', '--output', help='output_prefix', required=True)

    args = parser.parse_args()

    # parse the pan-genome file
    pangenes = parse_pangenes(args.pangenes)
    gff = parse_gff(args.gff)

    # add the synteny category
    pangenes_pivot = add_synteny_category(pangenes)
    # merge the pangenes and gff files
    gff_pangenes = merge_pangenes_gff(pangenes_pivot, gff)
    # select the syntelogs
    
    # make a pie chart
    make_pie_chart(gff_pangenes, args.output)

    # in the gff dataframe please replace the atttribute column for the rows that contain mRNA with the attribute column from the gff_pangenes dataframe, all other rows should remain the same, and the row from the gff dataframe and the gff should match on all columns except the attribute column. If the attribute column is empty in the gff_pangenes dataframe, the attribute column in the gff dataframe should kept as is.

    # add the synteny category and syntenic_genes to attributes collumn
    gff_pangenes['new_attributes'] = ';synteny_category=' + gff_pangenes['synteny_category'] + ';syntenic_genes=' + gff_pangenes['syntenic_genes']
    # strip tabs from the attributes column
    gff_pangenes['new_attributes'] = gff_pangenes['new_attributes'].str.replace('\t', '')

    # drop the columns that are not needed
    gff_pangenes.drop(columns=['synteny_category', 'syntenic_genes', 'haplotype', 'Synt_id'], inplace=True)
  
    # merge the two dataframes based on ID
    gff = pd.merge(gff, gff_pangenes, on=['ID', 'seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase'], how='left')
 
    
    # drop ID column
    print(gff['attributes_x'])
    print(gff['attributes_y'])

    gff['attributes'] = gff.apply(combine_attributes, axis=1)
    gff['attributes'] = gff['attributes'].str.replace('\t', '')

    gff.drop(columns=['ID', 'attributes_x', 'attributes_y', 'new_attributes'], inplace=True)

    print(gff)

    # save the dataframe to a tsv file
    gff.to_csv(f'{args.output}_synteny.gff', sep='\t', index=False, header=False)

    # filter the 4x syntelogs
    syntelogs = filter_4x_syntelogs(pangenes_pivot)

    syntelogs_lengths_list = []

    for attribute in ['exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR']:
        # add the length to the syntelogs
        ref_lengths = get_transcript_lengths(args.gff, attribute)
        syntelogs_lengths = add_length_to_syntelogs(syntelogs, ref_lengths, attribute)
        syntelogs_lengths['haplotype_id'] = syntelogs_lengths['haplotype_id'].apply(lambda x: ','.join(map(str, x)) if isinstance(x, list) else str(x))
        syntelogs_lengths_list.append(syntelogs_lengths)
        # make a barplot
        make_barplot(syntelogs_lengths, attribute, args.output)



# merge the dataframes in the list based on the shared columns
syntelogs_lengths_all = pd.merge(syntelogs_lengths_list[0], syntelogs_lengths_list[1], on=['Synt_id', 'haplotype_id'], how='inner')
syntelogs_lengths_all = pd.merge(syntelogs_lengths_all, syntelogs_lengths_list[2], on=['Synt_id', 'haplotype_id'], how='inner')
# save the dataframe to a tsv file
syntelogs_lengths_all.to_csv(f'{args.output}_syntelogs_all.tsv', sep='\t', index=False)

# expand the haplotype_id column
syntelogs_lengths_all['haplotype_id'] = syntelogs_lengths_all['haplotype_id'].str.split(',')
syntelogs_lengths_all = syntelogs_lengths_all.explode('haplotype_id')

syntelogs_lengths_all.to_csv(f'{args.output}_syntelogs_all_exploded.tsv', sep='\t', index=False)
