import pandas as pd
import argparse
import matplotlib.pyplot as plt


def dict_to_tsv(this_dict):
    df = pd.DataFrame(columns = this_dict.keys(),  index= this_dict['1G'].keys())
    for hap, error_dict in this_dict.items():
        row = []
        for rows, values in error_dict.items():
            row.append(values)
        df[hap] = row
    print(df)
def mapping_statistics(mapping_file, annotation_file):
    statistics_dict = {}
    # Function to calculate the number of missmpaing transripts per haplotoye 
    total_genes = len(annotation_file['gene_id'])
    total_wrong = len(mapping_file['real transcript'])


    for haplotype in ['1G', '2G', '3G', '4G']:
        statistics_dict[haplotype] = {}
        total_hap = len(annotation_file[annotation_file['gene_id'].str.contains(haplotype)])
        statistics_dict[haplotype]['total genes'] = round(total_hap / total_genes , 5)
        filtered_mapping = mapping_file[mapping_file['real transcript'].str.contains(haplotype)]
        wrong_mapping_total = len(mapping_file[mapping_file['real transcript'].str.contains(haplotype)])
        statistics_dict[haplotype]['wrong_mapping_total_%_' ] = round(wrong_mapping_total/total_wrong, 5)
        for haplotype2 in ['1G', '2G', '3G', '4G']:
            statistics_dict[haplotype]['wrong_mapping_to' + haplotype2] = round((len(filtered_mapping[filtered_mapping['mapped to'].str.contains(haplotype2)])) / total_wrong, 5)

    print(statistics_dict)
    return statistics_dict


# Set up command-line arguments
parser = argparse.ArgumentParser(description='Filter a CSV file to remove rows containing "," or "NA".')
parser.add_argument('mapping_file', help='TSV file that assignes reads to transcripts')
parser.add_argument('annotation_file', help='TSV file with gene location')


args = parser.parse_args()
mapping_file = args.mapping_file
annotation_file = args.annotation_file


mapping = pd.read_csv(mapping_file, sep = " ", header = None)
mapping.columns = ['real transcript', 'num', 'mapped to', 'num2']
mapping = mapping.drop(['num','num2'], axis=1)

annotation = pd.read_csv(annotation_file, sep = " ", header = None)

annotation.columns = ['gene_id', 'chrom', 'start', 'stop']

my_dict = mapping_statistics(mapping, annotation)
dict_to_tsv(my_dict)







