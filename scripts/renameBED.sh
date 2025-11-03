#!/bin/bash

# Check if an input file is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_file> <output_file>"
    exit 1
fi

# Assign the first argument to input
input="$1"
output="$2"

# Define two arrays
hap_numbers=("1" "2" "3" "4")
hap_letters=("a" "b" "c" "d")

prefix="a"

# Loop through the indices of the arrays
for i in "${!hap_numbers[@]}"; do
    # Get the elements from each array
    number="${hap_numbers[$i]}"
    letter="${hap_letters[$i]}"

    # Process the input file
    cat "$input" | grep -E "chr.*_$number" | \
    sed "s/chr_\([0-9]*\)_$number/$prefix$letter\1/g" | \
    sed "s/${letter}0/$letter/g" >> $output.tmp
done


# change the order of the bed file columns to match the format of the bed file used in the pipeline
awk '{print $1, $4, $2, $3}' OFS="\t" $output.tmp > $output
