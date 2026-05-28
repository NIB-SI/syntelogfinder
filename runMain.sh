
# Running the main pipeline with Nextflow, using the specified parameters and configuration.

nextflow run main.nf -resume -params-file params/my_example.json -c conf/nextflow.config -profile singularity --run_blast --mcscanx_path /MCScanX