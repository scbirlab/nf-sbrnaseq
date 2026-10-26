#!/usr/bin/env bash
# Usage: bash test/run-tests-illumina.sh [gh|crick|local] [ecoli|barnyard]
set -euox pipefail
HERE="$(readlink -f "$(dirname "${BASH_SOURCE[0]}")")"
REPO="$(readlink -f  "$HERE/..")"
PROFILE="${1:-gh}"
EXPT="${2:-ecoli}"
MM="${3:-no}"

case "$PROFILE" in
  gh)    
    export NXF_CONTAINER_ENGINE=docker
    EXTRA="-with-docker"
    ;;
  crick) 
    export NXF_CONTAINER_ENGINE=singularity; 
    PROFILE=local
    EXTRA="-with-singularity --eggnogdb ""$HERE""/_db/eggnogdb"
    module load Singularity Nextflow
    ;;
  local) EXTRA="-with-conda" ;;
  *) echo 'Usage: run-tests-illumina.sh [gh|crick|local] [ecoli|barnyard]' >&2; exit 2 ;;
esac

case "$EXPT" in
  ecoli) : ;;
  bsub) : ;;
  barnyard) : ;;
  *) echo 'Experiment must be ecoli, bsub, or barnyard' >&2; exit 2 ;;
esac

if [ "$EXPT" == "multimap" ]
then
  EXPT="$EXPT"-multimap
  MM_FLAG="--allow_multimap"
else
  MM_FLAG=
fi
mkdir -p "$HERE/illumina/$EXPT"
sheet="$HERE/illumina/$EXPT/sample-sheet.csv"
outdir="$HERE/illumina/$EXPT/outputs"

if [ ! -d "$HERE/illumina/$EXPT/fastq/" ]
then
  GENERATOR=(python3 "$HERE/scripts/generate_fixture.py" --mode "$EXPT" --output-dir "$HERE/illumina/$EXPT")
  if [[ -n "${SBRNASEQ_REFERENCE_DIR:-}" ]]; then
    GENERATOR+=(--reference-dir "$SBRNASEQ_REFERENCE_DIR")
  fi
  if [[ "${SBRNASEQ_NO_DOWNLOAD:-0}" == 1 ]]; then
    GENERATOR+=(--no-download)
  fi
  "${GENERATOR[@]}"
fi

nextflow run "$REPO" \
    -resume \
    -profile "$PROFILE" $EXTRA \
    -work-dir "$HERE/illumina/$EXPT/work" \
    -with-dag "$HERE/illumina/$EXPT/dag.html" \
    --mapper star "$MM_FLAG"\
    --sample_sheet "$sheet" \
    --inputs "$HERE/illumina/$EXPT/inputs" \
    --fastq_dir "$HERE/illumina/$EXPT/fastq" \
    --outputs "$outdir" \
    --trim_qual 5 \
    --min_length '9:38' \
    --end5_enabled \
    --end5_read 1 \
    --end5_boundary 5 \
    --min_counts_per_cell 10 \
    --strand 1
