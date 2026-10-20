
set -euox pipefail

HERE="$(dirname "$0")"
OUT="$HERE"/../illumina/barnyard/outputs/modifications
OUT=$(readlink -f "$OUT")
FILE=("$OUT"/*.mod-sites.tsv.gz)

for f in "${FILE[@]}"
do
    test -n "$f"

    python - "$f" <<'EOF'
import gzip
import csv
import sys

with gzip.open(sys.argv[1], "rt") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

assert rows, "No covered cell-site rows"

for r in rows:
    depth = int(r["callable_depth"])
    bases = sum(int(r[f"count_{b}"]) for b in "ACGT")
    assert depth == bases

    assert (
        int(r["ref_count"])
        + int(r["alt_count"])
        + int(r["other_count"])
        == depth
    )

    primary = int(r["mapping_primary_reads"])
    confident = int(r["mapping_confident_reads"])
    assert 0 <= confident <= primary

    if primary:
        assert abs(
            float(r["mapping_confident_fraction"])
            - confident / primary
        ) < 1e-8

print(f"PASS: {sys.argv[1]} -- {len(rows)} cell-site rows")
EOF

done
