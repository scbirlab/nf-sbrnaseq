
#!/usr/bin/env python
"""Regression tests for modification-site AnnData and MuData."""

from pathlib import Path
import subprocess
import sys
import tempfile

import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from test_utils.files import find_repo_root


BUILD_MOD = find_repo_root(Path(__file__), "bin/anndata/build-mod-anndata.py")
BUILD_MU = find_repo_root(Path(__file__), "bin/anndata/build-mudata.py")


def run(script, *args, success=True):
    result = subprocess.run(
        [sys.executable, str(script), *map(str, args)],
        capture_output=True,
        text=True,
    )
    if success and result.returncode:
        raise AssertionError(result.stdout + result.stderr)
    if not success and result.returncode == 0:
        raise AssertionError("Invalid input was accepted")
    return result


def dense(matrix):
    return matrix.toarray() if hasattr(matrix, "toarray") else np.asarray(matrix)


rows = [
    # Covered REF only.
    ("A", "site58", 5, 0, 0, 0),
    # Mixed REF/ALT.
    ("B", "site58", 3, 0, 2, 0),
    # No coverage, even though the row exists.
    ("C", "site58", 0, 0, 0, 0),
    # OTHER-only calls: covered, but ref_alt_depth=0.
    ("B", "site37", 0, 4, 0, 0),
]

metadata = {
    "site58": ("gene_trna1", "A", "G", "m1A58"),
    "site37": ("gene_trna2", "G", "A", "m1G37"),
}

records = []
for cell, site, a, c, g, t in rows:
    feature, ref_base, alt_base, modification = metadata[site]
    bases = dict(zip("ACGT", (a, c, g, t)))
    ref = bases[ref_base]
    alt = bases[alt_base]
    depth = sum(bases.values())

    records.append({
        "sample_id": "sample1",
        "cell_barcode": cell,
        "site_id": site,
        "feature_id": feature,
        "ref_base": ref_base,
        "alt_base": alt_base,
        "modification": modification,
        **{f"count_{b}": bases[b] for b in "ACGT"},
        "ref_count": ref,
        "alt_count": alt,
        "other_count": depth - ref - alt,
        "callable_depth": depth,
        "alt_fraction": (
            alt / (ref + alt) if ref + alt else np.nan
        ),
    })

df = pd.DataFrame(records)

with tempfile.TemporaryDirectory() as directory:
    tmp = Path(directory)
    counts_file = tmp / "counts.tsv"
    mod_file = tmp / "mod.h5ad"

    df.to_csv(counts_file, sep="\t", index=False)

    run(
        BUILD_MOD,
        "--counts", counts_file,
        "--sample-id", "sample1",
        "--output", mod_file,
    )

    mod = ad.read_h5ad(mod_file)
    assert mod.X.format == "csr"
    assert mod.var.loc["site58", "modification"] == "m1A58"
    assert mod.var.loc["site37", "modification"] == "m1G37"

    for base in "ACGT":
        assert f"count_{base}" in mod.layers

    base_total = sum(
        mod.layers[f"count_{b}"] for b in "ACGT"
    )
    np.testing.assert_array_equal(
        dense(base_total),
        dense(mod.layers["callable_depth"]),
    )

    combined = (
        mod.layers["ref_count"]
        + mod.layers["alt_count"]
        + mod.layers["other_count"]
    )
    np.testing.assert_array_equal(
        dense(combined),
        dense(mod.layers["callable_depth"]),
    )

    np.testing.assert_array_equal(
        dense(mod.X),
        dense(mod.layers["alt_count"]),
    )

    def value(layer, cell, site):
        i = mod.obs_names.get_loc(cell)
        j = mod.var_names.get_loc(site)
        return int(mod.layers[layer][i, j])

    assert value("callable_depth", "C", "site58") == 0
    assert value("alt_count", "C", "site58") == 0
    assert value("ref_count", "A", "site58") == 5
    assert value("alt_count", "A", "site58") == 0
    assert value("alt_count", "B", "site58") == 2
    assert value("ref_alt_depth", "B", "site58") == 5
    assert value("callable_depth", "B", "site37") == 4
    assert value("ref_alt_depth", "B", "site37") == 0

    print("PASS: base and REF/ALT count conservation")
    print("PASS: coverage and missingness semantics")
    print("PASS: modification metadata")

    # Corrupt one base count without updating depth.
    invalid = df.copy()
    invalid.loc[0, "count_A"] += 1
    bad_file = tmp / "invalid.tsv"
    invalid.to_csv(bad_file, sep="\t", index=False)

    result = run(
        BUILD_MOD,
        "--counts", bad_file,
        "--sample-id", "sample1",
        "--output", tmp / "invalid.h5ad",
        success=False,
    )
    assert "Base-count conservation failed" in result.stderr

    # Corrupt REF/ALT partition while preserving A/C/G/T.
    invalid = df.copy()
    invalid.loc[0, "ref_count"] -= 1
    invalid.to_csv(bad_file, sep="\t", index=False)

    result = run(
        BUILD_MOD,
        "--counts", bad_file,
        "--sample-id", "sample1",
        "--output", tmp / "invalid.h5ad",
        success=False,
    )
    assert "REF/ALT count conservation failed" in result.stderr

    print("PASS: invalid counts rejected")

    # Build RNA and filtered RNA with partially overlapping cells.
    rna = ad.AnnData(
        X=csr_matrix([[10, 2], [4, 6], [1, 3]]),
        obs=pd.DataFrame(index=["A", "B", "D"]),
        var=pd.DataFrame(index=["gene1", "gene2"]),
    )
    rna.layers["counts"] = rna.X.copy()

    rna_file = tmp / "rna.h5ad"
    qc_file = tmp / "rna-qc.h5ad"
    rna.write_h5ad(rna_file)
    rna[["A", "B"], :].copy().write_h5ad(qc_file)

    output = tmp / "sample.h5mu"
    run(
        BUILD_MU,
        "--rna", rna_file,
        "--rna-qc", qc_file,
        "--modalities", mod_file,
        "--sample-id", "sample1",
        "--output", output,
    )

    mdata = mu.read_h5mu(output)
    joined = mdata.mod["mod"]

    assert list(mdata.obs_names) == [
        "sample1:A", "sample1:B",
        "sample1:C", "sample1:D",
    ]
    assert joined.var.loc["site58", "modification"] == "m1A58"

    for layer in mod.layers.keys():
        assert int(joined.layers[layer].sum()) == int(
            mod.layers[layer].sum()
        ), layer

    assert mdata.obs.loc["sample1:C", "has_mod"]
    assert not mdata.obs.loc["sample1:D", "has_mod"]
    assert not mdata.obs.loc["sample1:D", "rna_qc_pass"]

    # D is structurally absent from modification measurements.
    assert joined["sample1:D", :].X.sum() == 0
    assert not joined.obs.loc[
        "sample1:D", "present_in_modality"
    ]

    print("PASS: MuData count and metadata preservation")
    print("PASS: structural absence and barcode alignment")
