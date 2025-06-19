# potato-allelic-orthogroups

Nextflow pipeline to group genes on polyploid phased assemblies that are orthologous and syntelogous based on GENESPACE.

minimal input:
- genome fasta of phased reference
- gff with CDS corresponding to the reference

config file should look like this


{
    "reference_fasta": "genome.fa",
    "reference_gff": "annotation.gff3",
    "outdir": "output_path"
}




## Output

 category grouping




