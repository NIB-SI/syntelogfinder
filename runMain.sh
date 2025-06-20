
# Desiree with desiree annotation

nextflow run main.nf -resume -params-file params/desiree_liftoff.json -c conf/nextflow.config -profile conda          



nextflow run main.nf -resume -params-file params/desiree.json -with-conda -c conf/nextflow.config -profile conda   


nextflow run main.nf -resume -params-file params/atlantic.json -with-conda -c conf/nextflow.config -profile conda  