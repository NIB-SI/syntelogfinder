# potato-allelic-orthogroups




anaylsing transcriptome complexitiy:
- download http://spuddb.uga.edu/data/ATL_v3/ATL_v3.hc_gene_models.gff3.gz
- gunzip and move to gff_references


Nextflow pipeline to unify gene IDs of phased potato references based on orthogroups 

pipeline to create input for SAynVisio to vizualize the 12 x 4n chromosomes from Otava, Atlantic and Castle Russet

Codong sequence or transcript sequence?


Input_Proteoms

Otava downloaded here:
http://spuddb.uga.edu/data/St2.prot.fa.gz
http://spuddb.uga.edu/data/St1.prot.fa.gz
http://spuddb.uga.edu/data/He1.prot.fa.gz
http://spuddb.uga.edu/data/He2.prot.fa.gz
cat * > Otava.fa


Atlantic
high confidence gene models





http://spuddb.uga.edu/data/ATL_v3/ATL_v3.working_models.repr.pep.fa.gz
 mv ATL_v3.working_models.repr.pep.fa.gz atlantic.fa


castle russet
http://spuddb.uga.edu/data/cr.hc.repr.pm.pep.fasta.gz



Downlad genome fasta and gff to potato_culitvar_references for 
- atlantic


Atalntic chrom names:

>chr01_1 chrom 1 hap 1
>chr01_2
>chr01_3
>chr01_4 
>chr01_0 chrom 1 hap unnasigned?
>chr02_1
>chr02_2
>chr02_3

example gene name on chr11_1: Soltu.Atl_v3.11_1G019420.1

Otava
He1-St01

split into haplotypes
(nextflow) egret:nadjafn$ seqkit grep -r -p 'Soltu.Atl_v3.*_2' atlantic.fa > atlantic_2.fa
(nextflow) egret:nadjafn$ seqkit grep -r -p 'Soltu.Atl_v3.*_1' atlantic.fa > atlantic_1.fa
(nextflow) egret:nadjafn$ seqkit grep -r -p 'Soltu.Atl_v3.*_4' atlantic.fa > atlantic_4.fa

(nextflow) egret:nadjafn$ seqkit grep -r -p 'Soltu.Cru.*_2' castle_russet.fa > castle_russet_2.fa
(nextflow) egret:nadjafn$ seqkit grep -r -p 'Soltu.Cru.*_3' castle_russet.fa > castle_russet_3.fa
(nextflow) egret:nadjafn$ seqkit grep -r -p 'Soltu.Cru.*_1' castle_russet.fa > castle_russet_1.fa
(nextflow) egret:nadjafn$ seqkit grep -r -p 'Soltu.Cru.*_4' castle_russet.fa > castle_russet_4.fa


downlaoad gff files
Modify chrom names for scScanX
grep '^Chr' He1.protein-coding.gene.gff3 | sed 's#^#Otava_He1#' > Otava_He1.gff3 
grep '^chr' CR_v2.0.hc_gene_models.repr.gff3 |  sed 's#^# Castle_russet_#' > Castle_russet.gff3
grep '^chr' ATL_v3.hc_gene_models.repr.gff3 |  sed 's#^#Atlantic_#' > Atlantic.gff3 

cat Atlantic.gff3 Castle_russet.gff3 Otava_*.gff3 > cultivars_merged.gff3
gff2bed --do-not-sort  < cultivars_merged.gff3 > cultivars_merged.bed

/scratch/nadjafn/MCScanX/MCScanX Input_MCscanX/cultivars

--> no scaffolds and hap 0 (unassigned?)


check 0e/f0c4fa

Processes neccesary:


Before Running: 


download and put into Input:

Castle russet repr protem [http://spuddb.uga.edu/data/cr.hc.repr.pm.pep.fasta.gz]

Atlantic pepr protem []



1) downladed prot and gff files
3) split prot into haplotpyes --> run orthofinder
2) merge all prot create diamond db --> run diamond
3) change gff file chrom name to culitvar_chromname, merge gff3, convert to bed
4) for mcxscan put diamond as culitvar.blast and cultivar.gff in same dir
5) run_MCScan

