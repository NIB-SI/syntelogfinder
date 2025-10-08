
# Desiree with desiree annotation

nextflow run main.nf -resume -params-file params/desiree_liftoff.json -c conf/nextflow.config -profile conda
nextflow run main.nf -resume -params-file params/desiree_liftoff_triploid.json -c conf/nextflow.config -profile conda


nextflow run main.nf -resume -params-file params/desiree.json -with-conda -c conf/nextflow.config -profile conda


nextflow run main.nf -resume -params-file params/atlantic.json -with-conda -c conf/nextflow.config -profile conda

# Atlantic with liftoff
nextflow run main.nf -resume -params-file params/atlantic_liftoff.json -c conf/nextflow.config -profile conda


# Desiree with liftoff and Bambu
nextflow run main.nf -resume -params-file params/desiree_liftoff_bambu.json -c conf/nextflow.config -profile conda


# SWEEEET
nextflow run main.nf -resume -params-file params/sweetpotato.json -c conf/nextflow.config -profile conda

# wheat
nextflow run main.nf -resume -params-file params/wheatAK58.json -c conf/nextflow.config -profile conda --run_blast


# rice
nextflow run main.nf -resume -params-file params/rice_Nip.json -c conf/nextflow.config -profile conda --run_blast --mcscanx_path /DKED/scratch/nadjafn/MCScanX
