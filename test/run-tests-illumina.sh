 #!/usr/bin/env bash

set -euox pipefail

script_dir="$(dirname $0)"
DOCKER=${1:-no}

if [ "$DOCKER" == "gh" ]
then
    export NXF_CONTAINER_ENGINE=docker
    docker_flag='-profile gh'
elif [ "$DOCKER" == "crick" ]
then
    export NXF_CONTAINER_ENGINE=singularity
    module load Nextflow Singularity
    docker_flag='-with-singularity -profile local'
else
    export SINGULARITY_FAKEROOT=1
    docker_flag=''
fi

nextflow run "$script_dir"/.. \
    -resume $docker_flag \
    -work-dir "$script_dir"/illumina/work \
    -with-dag "$script_dir"/illumina/dag.html \
    --mapper star \
    --sample_sheet "$script_dir"/illumina/sample-sheet.csv \
    --inputs "$script_dir"/illumina/inputs \
    --fastq_dir "$script_dir"/illumina/fastq \
    --outputs "$script_dir"/illumina/outputs
