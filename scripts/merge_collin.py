import pandas as pd
import argparse

def parse_datagram(lines):
    # Extracting the header information
    header = lines[0].split()
    alignment_id = header[2].strip(':')
    score = float(header[3].split('=')[1])
    e_value = header[4].split('=')[1]
    n = int(header[5].split('=')[1])
    
    # Parsing each line to extract the data
    data = []
    for line in lines[1:]:
        parts = line.split()
        index = parts[0].split('-')[1]
        gene1 = parts[2]
        gene2 = parts[3]
        e_value_line = parts[3] if len(parts) == 4 else None

        # Check if either of the gene IDs contains 'Soltu.Atl_v3.'
        #if 'Soltu.Atl_v3.' in gene1 and 'Soltu.Atl_v3.' in gene2:
        data.append([alignment_id, index, gene1, gene2, e_value_line])
    
    # Creating a DataFrame
    if data:  # Only create a DataFrame if there is data
        df = pd.DataFrame(data, columns=['Alignment', 'Index', 'Gene1', 'Gene2', 'E_Value_Line'])
        df['Score'] = score
        df['E_Value'] = e_value
        df['N'] = n
        return df
    else:
        return None

# Setup the argument parser
parser = argparse.ArgumentParser(description='Process datagrams from a file.')
parser.add_argument('file_path', help='Path to the data file')
args = parser.parse_args()

# Reading the file and splitting into datagrams
with open(args.file_path, 'r') as file:
    lines = file.readlines()

datagrams = []
current_datagram = []

# skip the lines until first appearance of '## '
while not lines[0].startswith('## '):
    lines.pop(0)

for line in lines:
    if line.startswith('## '):
        if current_datagram:
            df = parse_datagram(current_datagram)
            if df is not None:
                datagrams.append(df)
            current_datagram = []
        current_datagram.append(line.strip())
    elif line.strip():
        current_datagram.append(line.strip())

# Add the last datagram if it exists
if current_datagram:
    df = parse_datagram(current_datagram)
    if df is not None:
        datagrams.append(df)

# Concatenate all DataFrames
final_df = pd.concat(datagrams, ignore_index=True) if datagrams else pd.DataFrame()

# Display the final DataFrame
print(final_df)
