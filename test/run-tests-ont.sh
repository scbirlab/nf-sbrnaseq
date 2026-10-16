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
    -work-dir "$script_dir"/ont/work \
    -with-dag "$script_dir"/ont/dag.html \
    --nanopore \
    --mapper minimap2 \
    --reverse '2 3' \
    --sample_sheet "$script_dir"/ont/sample-sheet.csv \
    --inputs "$script_dir"/ont/inputs \
    --fastq_dir "$script_dir"/ont/fastq \
    --outputs "$script_dir"/ont/outputs
