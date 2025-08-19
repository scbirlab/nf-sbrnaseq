 #!/usr/bin/env bash

set -x
set -e

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

cd "$script_dir"/ont
nextflow run "$script_dir"/../.. \
    -resume $docker_flag \
    --sample_sheet "$script_dir"/ont/sample-sheet.csv \
    --inputs "$script_dir"/ont/inputs \
    --outputs "$script_dir"/ont/outputs
cd ..
