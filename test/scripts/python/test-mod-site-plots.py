
#!/usr/bin/env python3
"""Validate site-specific modification plotting outputs."""

import subprocess
import sys
import tempfile
from pathlib import Path

import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from test_utils.files import find_repo_root


SCRIPT = find_repo_root(Path(__file__), "bin/viz/plot-mod-sites.py")

with tempfile.TemporaryDirectory() as directory:
    tmp = Path(directory)

    barcodes = ["A", "B", "C", "D"]

    rna = ad.AnnData(
        X=csr_matrix(np.ones((4, 3))),
        obs=pd.DataFrame({
            "sample_id": ["sample1"] * 4,
            "leiden": pd.Categorical(["0", "0", "1", "1"]),
        }, index=barcodes),
    )
    rna.obsm["X_umap"] = np.array([
        [0., 0.], [1., 0.],
        [0., 1.], [1., 1.],
    ])
    rna_path = tmp / "cluster.h5ad"
    rna.write_h5ad(rna_path)

    mod = ad.AnnData(
        X=csr_matrix([[0], [0], [2], [0]]),
        obs=pd.DataFrame({
            "present_in_modality": [True, True, True, False],
        }, index=[f"sample1:{b}" for b in barcodes]),
        var=pd.DataFrame(index=["site58"]),
    )
    mod.layers["ref_count"] = csr_matrix(
        [[0], [5], [3], [0]]
    )
    mod.layers["alt_count"] = csr_matrix(
        [[0], [0], [2], [0]]
    )
    mod.layers["callable_depth"] = csr_matrix(
        [[0], [5], [5], [0]]
    )
    mod.layers["ref_alt_depth"] = csr_matrix(
        [[0], [5], [5], [0]]
    )

    mdata = mu.MuData({"mod": mod})
    mdata_path = tmp / "sample.h5mu"
    mdata.write_h5mu(mdata_path)

    subprocess.run([
        sys.executable, str(SCRIPT),
        "--clustered", str(rna_path),
        "--mudata", str(mdata_path),
        "--output", str(tmp / "modsites"),
    ], check=True)

    coverage = pd.read_csv(
        tmp / "modsites.cluster-coverage.csv",
        index_col=0,
    )
    fractions = pd.read_csv(
        tmp / "modsites.cluster-alt-fraction.csv",
        index_col=0,
    )

    assert np.isclose(coverage.loc["site58", "0"], 0.5)
    assert np.isclose(coverage.loc["site58", "1"], 1.0)
    assert np.isclose(fractions.loc["site58", "0"], 0.0)
    assert np.isclose(fractions.loc["site58", "1"], 0.4)

    assert (
        tmp / "modsites.umap-coverage.page01.png"
    ).exists()
    assert (
        tmp / "modsites.umap-fraction.page01.png"
    ).exists()

    print("PASS: coverage and ALT fractions")
    print("PASS: uncovered versus REF-only cells")
    print("PASS: modification UMAP figures")
