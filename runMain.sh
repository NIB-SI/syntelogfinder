
# Desiree with desiree annotation

nextflow run main.nf -resume -params-file params/desiree_liftoff.json -c conf/nextflow.config -profile conda --run_blast --mcscanx_path /DKED/scratch/nadjafn/MCScanX
nextflow run main.nf -resume -params-file params/my_example.json -c conf/nextflow.config -profile conda --run_blast --mcscanx_path /DKED/scratch/nadjafn/MCScanX


nextflow run main.nf -resume -params-file params/desiree.json -with-conda -c conf/nextflow.config -profile conda


nextflow run main.nf -resume -params-file params/atlantic.json -with-conda -c conf/nextflow.config -profile conda

# Atlantic with liftoff
nextflow run main.nf -resume -params-file params/atlantic_with_liftoff.json -c nextflow.config -profile conda --run_blast --mcscanx_path /DKED/scratch/nadjafn/MCScanX



# Desiree with liftoff and Bambu
nextflow run main.nf -resume -params-file params/desiree_liftoff_bambu.json -c conf/nextflow.config -profile conda --run_blast --mcscanx_path /DKED/scratch/nadjafn/MCScanX


# SWEEEET
nextflow run main.nf -resume -params-file params/sweetpotato.json -c conf/nextflow.config -profile conda

# wheat
nextflow run main.nf -resume -params-file params/wheatAK58.json -c conf/nextflow.config -profile conda --run_blast


# rice
nextflow run main.nf -resume -params-file params/rice_Nip.json  -profile singularity --run_blast --mcscanx_path /DKED/scratch/nadjafn/MCScanX




nextflow run main.nf -resume -params-file params/my_example.json -c conf/nextflow.config -profile conda --run_blast --mcscanx_path /DKED/scratch/nadjafn/MCScanX