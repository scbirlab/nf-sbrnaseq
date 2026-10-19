#!/usr/bin/env bash
# Usage: bash test/run-tests-illumina.sh [gh|crick|local] [ecoli|barnyard]
set -euox pipefail
HERE="$(readlink -f "$(dirname "${BASH_SOURCE[0]}")")"
REPO="$(readlink -f  "$HERE/..")"
PROFILE="${1:-gh}"
EXPT="${2:-ecoli}"

case "$PROFILE" in
  gh)    
    export NXF_CONTAINER_ENGINE=docker
    EXTRA="-with-docker"
    ;;
  crick) 
    export NXF_CONTAINER_ENGINE=singularity; 
    PROFILE=local
    EXTRA="-with-singularity"
    module load Singularity Nextflow
    ;;
  local) EXTRA="-with-conda" ;;
  *) echo 'Usage: run-tests-illumina.sh [gh|crick|local] [ecoli|barnyard]' >&2; exit 2 ;;
esac
sheet="$HERE/illumina/$EXPT/sample-sheet.csv"
outdir="$HERE/illumina/$EXPT/outputs"
case "$EXPT" in
  ecoli) : ;;
  bsub) : ;;
  barnyard) : ;;
  *) echo 'Experiment must be ecoli, bsub, or barnyard' >&2; exit 2 ;;
esac

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
    -profile "$PROFILE" "$EXTRA" \
    -work-dir "$HERE/illumina/$EXPT/work" \
    -with-dag "$HERE/illumina/$EXPT/dag.html" \
    --mapper star \
    --sample_sheet "$sheet" \
    --inputs "$HERE/illumina/$EXPT/inputs" \
    --fastq_dir "$HERE/illumina/$EXPT/fastq" \
    --outputs "$outdir" \
    --trim_qual 5 \
    --min_length '9:38' \
    --mod_sites "$HERE/illumina/trna-mod/trna-sites.tsv" \
    --sites_join_on feature \
    --strand 1
