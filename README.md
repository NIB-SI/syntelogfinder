# potato-allelic-orthogroups

Nextflow pipeline to group genes on polyploid phased assemblies that are orthologous and syntelogous based on GENESPACE.

minimal input:
- genome fasta of phased reference
- gff with CDS corresponding to the reference

config file should look like this


{
    "reference_fasta": "genome.fa",
    "reference_gff": "annotation.gff",
    "outdir": "output_path"
}

Run like this:
nextflow run main.nf -resume -params-file params/params.json -c cond/nextflow.config -profile conda


the gff file should look like this https://agat.readthedocs.io/en/latest/gff_to_gtf.html#the-gff-file-to-convert

gene, mRNA, exon, CDS features



## Output

 category grouping




