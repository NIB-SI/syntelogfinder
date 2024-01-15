import pandas as pd
import argparse

def filter_csv(input_file, output_file):
    # Read the CSV file
    data = pd.read_csv(input_file, sep = "\t")
    print(data)
    filtered_data = data[~data.astype(str).apply(lambda row: row.str.contains(',')).any(axis=1)]
    filtered_data = filtered_data.dropna(how='any')
    print(filtered_data)
    # Filter Orthogroups that dont have gene on each hpalotype

 
    # # Create a condition to filter rows
    # condition = data.apply(lambda row: row.str.startswith('/s')).sum(axis=1) == 4
    # # Extract rows and select specific columns
    # result_df = data.loc[condition, data.columns[data.columns.str.startswith('Soltu.Atl') | (data.columns == 'Orthogroup')]]

   
    # print(result_df)
    # # Filter rows that don't contain ',' or 'NA' in any column
    # filtered_data = data[~data.astype(str).apply(lambda row: row.str.contains(',')).any(axis=1)]
    # # Filter rows without any value in any column
    # filtered_data = filtered_data.dropna(how='any')
    
    # # Write the filtered data to a new CSV file
    # filtered_data.to_csv(output_file, index=False,sep = "\t")
    # print(len(filtered_data))

# Set up command-line arguments
parser = argparse.ArgumentParser(description='Filter a CSV file to remove rows containing "," or "NA".')
parser.add_argument('input_file', help='Input CSV file')
parser.add_argument('output_file', help='Output filtered CSV file')

args = parser.parse_args()
input_file = args.input_file
output_file = args.output_file

# Filter the CSV file
filter_csv(input_file, output_file)