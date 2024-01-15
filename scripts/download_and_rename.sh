#!/bin/bash

# Function to download and check for errors
download() {
    url=$1
    output=$2
    wget -O "$output" "$url" || { echo "Failed to download $url"; exit 1; }
}

# Create directories
mkdir -p Haplotype_Split GFF_Split || { echo "Failed to create directories"; exit 1; }

# Process CR
mkdir -p Input/cr
(
    cd Input/cr
    download http://spuddb.uga.edu/data/cr.hc.repr.pm.pep.fasta.gz "cr.hc.repr.pm.pep.fasta.gz"
    download http://spuddb.uga.edu/data/cr.hc.repr.pm.locus_assign.gff3.gz "cr.hc.repr.pm.locus_assign.gff3.gz"

        # Define two arrays
    hap_numbers=("1" "2" "3" "4")
    hap_letters=("a" "b" "c" "d")
    
    # Check if both arrays have the same length
    if [ "${#hap_numbers[@]}" -ne "${#hap_letters[@]}" ]; then
        echo "Arrays do not have the same length."
        exit 1
    fi

    # Loop through the indices of the arrays
    for i in "${!hap_numbers[@]}"; do
        # Get the elements from each array
        number="${hap_numbers[$i]}"
        letter="${hap_letters[$i]}"

        zcat cr.hc.repr.pm.pep.fasta.gz | seqkit grep -r -p "Soltu.Cru.*_$number" > "../../Haplotype_Split/C$letter.fa"
        zcat cr.hc.repr.pm.locus_assign.gff3.gz | grep -E "chr.*_$number" | sed  "s/chr\([0-9]*\)_$number/C$letter\1/g" | sed "s/C${letter}0/C$letter/g" | awk -F'\t' '$1 !~ /UA/' > "../../GFF_Split/C$letter.gff"
    done
)

# Process ATL
mkdir -p Input/at
(
    cd Input/at
    download http://spuddb.uga.edu/data/ATL_v3/ATL_v3.hc_gene_models.repr.pep.fa.gz "ATL_v3.hc_gene_models.repr.pep.fa.gz"
    download http://spuddb.uga.edu/data/ATL_v3/ATL_v3.hc_gene_models.repr.gff3.gz "ATL_v3.hc_gene_models.repr.gff3.gz"

     # Define two arrays
    hap_numbers=("1" "2" "3" "4")
    hap_letters=("a" "b" "c" "d")
    
    # Check if both arrays have the same length
    if [ "${#hap_numbers[@]}" -ne "${#hap_letters[@]}" ]; then
        echo "Arrays do not have the same length."
        exit 1
    fi

    # Loop through the indices of the arrays
    for i in "${!hap_numbers[@]}"; do
        # Get the elements from each array
        number="${hap_numbers[$i]}"
        letter="${hap_letters[$i]}"

        zcat ATL_v3.hc_gene_models.repr.pep.fa.gz | seqkit grep -r -p "Soltu.Atl_v3.*_$number" > "../../Haplotype_Split/A$letter.fa"
        zcat ATL_v3.hc_gene_models.repr.gff3.gz | grep -E "chr.*_$number" | sed  "s/chr\([0-9]*\)_$number/A$letter\1/g" | sed "s/A${letter}0/A$letter/g"  > "../../GFF_Split/A$letter.gff"
    done
)

# Process OT
mkdir -p Input/ot
(
    cd Input/ot
    for file in He1.prot He2.prot St1.prot St2.prot; do
        download "http://spuddb.uga.edu/data/$file.fa.gz" "$file.fa.gz"
    done
    zcat He1.prot.fa.gz > ../../Haplotype_Split/Oa.fa
    zcat He2.prot.fa.gz > ../../Haplotype_Split/Ob.fa
    zcat St1.prot.fa.gz > ../../Haplotype_Split/Oc.fa
    zcat St2.prot.fa.gz > ../../Haplotype_Split/Od.fa

    

    for file in He1.protein-coding.gene He2.protein-coding.gene St1.protein-coding.gene St2.protein-coding.gene; do
        download "http://spuddb.uga.edu/data/$file.gff3.gz" "$file.gff3.gz"
    done
    # rename chromosomes
    zcat He1.protein-coding.gene.gff3.gz | sed 's/Chr/Oa/g' | awk -F'\t' '$1 !~ /UA/' > ../../GFF_Split/Oa.gff
    zcat He2.protein-coding.gene.gff3.gz | sed 's/Chr/Ob/g' | awk -F'\t' '$1 !~ /UA/' > ../../GFF_Split/Ob.gff
    zcat St1.protein-coding.gene.gff3.gz | sed 's/Chr/Oc/g'| awk -F'\t' '$1 !~ /UA/' > ../../GFF_Split/Oc.gff
    zcat St2.protein-coding.gene.gff3.gz | sed 's/Chr/Od/g' | awk -F'\t' '$1 !~ /UA/' > ../../GFF_Split/Od.gff
)

echo "Processing complete."
