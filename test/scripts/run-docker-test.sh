
#!/usr/bin/env bash
set -euo pipefail

# Run nf-sbrnaseq tests inside Docker.
#
# Usage:
#   scripts/run-docker-test.sh test/scripts/python/test-mudata-alignment.py
#   scripts/run-docker-test.sh test/scripts/python/*.py
#   scripts/run-docker-test.sh -- python -c 'import scanpy; print(scanpy.__version__)'
#
# Override image:
#   IMAGE=ghcr.io/scbirlab/nf-sbrnaseq:latest scripts/run-docker-test.sh ...

REPO="$(readlink -f "$(dirname "${BASH_SOURCE[0]}")/../..")"
# IMAGE="${IMAGE:-nf-sbrnaseq:test}"
IMAGE="${IMAGE:-"ghcr.io/scbirlab/nf-sbrnaseq:latest"}"
PYTHON=/opt/conda/envs/env/bin/python

if [[ $# -eq 0 ]]; then
    echo "Usage: $0 TEST.py [TEST.py ...] | -- COMMAND [args ...]" >&2
    exit 1
fi

command -v docker >/dev/null 2>&1 || {
    echo "ERROR: Docker is not installed" >&2
    exit 1
}

# Same temporary filesystem semantics as the Singularity runner.
SCRATCH="$(mktemp -d)"
trap 'rm -rf "$SCRATCH"' EXIT

# Docker defaults to the image user (UID 1000).
# Use that user and permit writes to the temporary mount.
chmod 777 "$SCRATCH"

if [[ "$1" == "--" ]]; then
    shift

    if [[ $# -eq 0 ]]; then
        echo "ERROR: Missing command after --" >&2
        exit 1
    fi

    CMD=("$@")
else
    # Execute each Python test, stopping on the first failure.
    CMD=(
        /bin/bash -c
        'set -euo pipefail
         for test_file in "$@"; do
             echo "=== Running $test_file ==="
             /opt/conda/envs/env/bin/python "$test_file"
         done'
        bash
        "$@"
    )
fi

echo "Repository: $REPO"
echo "Image:      $IMAGE"
echo "Command:    ${CMD[*]}"
echo

docker run --rm \
    --user "$(id -u):$(id -g)" \
    -e MPLCONFIGDIR=/tmp/matplotlib-cache \
    -e NUMBA_CACHE_DIR=/tmp/numba-cache \
    --network none \
    --workdir /workspace \
    --mount "type=bind,source=$REPO,target=/workspace,readonly" \
    --mount "type=bind,source=$SCRATCH,target=/test-tmp" \
    --tmpfs /tmp:rw,mode=1777 \
    --env HOME=/test-tmp \
    --env TMPDIR=/test-tmp \
    --env XDG_CACHE_HOME=/test-tmp/xdg \
    --env MPLCONFIGDIR=/test-tmp/mpl \
    --env NUMBA_CACHE_DIR=/test-tmp/numba \
    --env MPLBACKEND=Agg \
    --env PATH=/opt/conda/envs/env/bin:/usr/local/bin:/usr/bin:/bin \
    --entrypoint /usr/bin/env \
    "$IMAGE" \
    "${CMD[@]}"
