
#!/usr/bin/env bash
set -euox pipefail

# Run tests or commands inside the nf-sbrnaseq Singularity image.
#
# Usage:
#   scripts/run-container-test.sh tests/test-mudata-alignment.py
#   scripts/run-container-test.sh tests/test-anndata-filtering.py
#   scripts/run-container-test.sh -- pytest tests/ -v
#   scripts/run-container-test.sh -- python -c 'import scanpy; print(scanpy.__version__)'
#
# Override image discovery:
#   IMAGE=/path/to/container.sif scripts/run-container-test.sh tests/test.py

REPO="$(readlink -f "$(dirname "${BASH_SOURCE[0]}")/../..")"
ENGINE="${ENGINE:-}"

if [[ -z "$ENGINE" ]]; then
    if command -v singularity >/dev/null 2>&1; then
        ENGINE=singularity
    elif command -v apptainer >/dev/null 2>&1; then
        ENGINE=apptainer
    else
        echo "ERROR: Neither singularity nor apptainer found" >&2
        exit 1
    fi
fi

if [[ $# -eq 0 ]]; then
    echo "Usage: $0 TEST.py [args...] | -- COMMAND [args...]" >&2
    exit 1
fi

# Discover cached image unless explicitly provided.
if [[ -z "${IMAGE:-}" ]]; then
    shopt -s nullglob
    images=(
        "$REPO"/.singularity/*nf-sbrnaseq*.sif
        "$REPO"/.singularity/*nf-sbrnaseq*.img
    )
    shopt -u nullglob

    if [[ ${#images[@]} -ne 1 ]]; then
        echo "ERROR: Expected exactly one cached nf-sbrnaseq image." >&2
        echo "Found ${#images[@]} matching images." >&2
        echo "Set IMAGE=/absolute/path/to/image.sif explicitly." >&2
        exit 1
    fi

    IMAGE="${images[0]}"
fi

# Resolve local image paths before entering the container.
if [[ "$IMAGE" != *://* ]]; then
    IMAGE="$(realpath "$IMAGE")"
    [[ -f "$IMAGE" ]] || {
        echo "ERROR: Image not found: $IMAGE" >&2
        exit 1
    }
fi

# Writable host-backed temporary directory.
SCRATCH="$(mktemp -d "${TMPDIR:-/tmp}/sbrnaseq-test.XXXXXXXX")"
trap 'rm -rf "$SCRATCH"' EXIT

# The Dockerfile installs Python inside this named environment.
PYTHON=/opt/conda/envs/env/bin/python
CONDA_BIN=/opt/conda/envs/env/bin

if [[ "$1" == "--" ]]; then
    shift
    [[ $# -gt 0 ]] || {
        echo "ERROR: Missing command after --" >&2
        exit 1
    }
    CMD=("$@")
else
    CMD=("$PYTHON" "$@")
fi

echo "Repository: $REPO"
echo "Image:      $IMAGE"
echo "Command:    ${CMD[*]}"
echo

"$ENGINE" exec \
    --cleanenv \
    --no-home \
    --bind "$REPO:/workspace" \
    --bind "$SCRATCH:/test-tmp" \
    --pwd /workspace \
    --env "PATH=$CONDA_BIN:/usr/local/bin:/usr/bin:/bin" \
    --env "TMPDIR=/test-tmp" \
    --env "XDG_CACHE_HOME=/test-tmp/xdg" \
    --env "MPLCONFIGDIR=/test-tmp/mpl" \
    --env "NUMBA_CACHE_DIR=/test-tmp/numba" \
    --env "MPLBACKEND=Agg" \
    "$IMAGE" \
    "${CMD[@]}"
