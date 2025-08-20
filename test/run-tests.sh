 #!/usr/bin/env bash

set -euox pipefail

script_dir="$(dirname $0)"
DOCKER=${1:-no}

if [ "$DOCKER" == "gh" ]
then
    export NXF_CONTAINER_ENGINE=docker
    docker_flag='-profile gh'
else
    export SINGULARITY_FAKEROOT=1
    docker_flag=''
fi

nextflow run "$script_dir"/.. \
    -resume $docker_flag \
    --nanopore \
    --mapper minimap2 \
    --reverse '2 3' \
    --sample_sheet "$script_dir"/ont/sample-sheet.csv \
    --inputs "$script_dir"/ont/inputs \
    --fastq_dir "$script_dir"/ont/fastq \
    --outputs "$script_dir"/ont/outputs
