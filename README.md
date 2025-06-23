# potato-allelic-orthogroups

Nextflow pipeline to group genes on polyploid phased assemblies that are orthologous and syntelogous based on GENESPACE.

Requirements:

- nextflow

- conda

The following packages are not in bioconda/pip so need to be installed manually:
- McxScan (follow instructions [here](/scratch/nadjafn/potato-allelic-orthogroups/modules/local/genespace/genespace_run/environment.yml) and provide path to installation to --mcscanx_path)
- GENESPACE ([instructions](https://github.com/jtlovell/GENESPACE?tab=readme-ov-file#2-software-installation))(in side the conda env genespace-env (potato-allelic-orthogroups/modules/local/genespace/genespace_run/environment.yml))

minimal input:
- genome fasta of phased reference
- gff with CDS corresponding to the reference

```
{
    "reference_fasta": "genome.fa",
    "reference_gff": "annotation.gff",
    "outdir": "output_path"
}
```

Run like this:
```
nextflow run main.nf -resume -params-file params/params.json -c cond/nextflow.config -profile conda --mcscanx_path [path to MCScanX folder]
```


the gff file should look like this https://agat.readthedocs.io/en/latest/gff_to_gtf.html#the-gff-file-to-convert

gene, mRNA, exon, CDS features



## Output

 category grouping




