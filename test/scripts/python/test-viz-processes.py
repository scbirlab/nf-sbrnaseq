
#!/usr/bin/env python
"""Smoke and numerical tests for independent exploratory analyses."""

from pathlib import Path
import subprocess
import sys
import tempfile

import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix



ROOT = next(
    p for p in Path(__file__).resolve().parents
    if (p / "bin/viz/embed-trna.py").exists()
)
BIN = ROOT / "bin/viz"


def run(script, *args):
    result = subprocess.run(
        [sys.executable, str(BIN / script), *map(str, args)],
        capture_output=True,
        text=True,
    )
    if result.returncode:
        raise AssertionError(result.stdout + result.stderr)


rng = np.random.default_rng(42)
n = 40
cells = [f"cell{i}" for i in range(n)]
qualified = [f"sample1:{c}" for c in cells]

signal = np.arange(1, n + 1)
counts = np.column_stack([
    signal,
    2 * signal,
    rng.integers(2, 30, (n, 5)),
])

rna = ad.AnnData(
    X=csr_matrix(np.log1p(counts)),
    obs=pd.DataFrame({
        "sample_id": "sample1",
        "leiden": pd.Categorical(["0"] * 20 + ["1"] * 20),
        "total_counts": counts.sum(axis=1),
        "pct_counts_is_tRNA": (
            counts[:, 2:].sum(axis=1) / counts.sum(axis=1) * 100
        ),
    }, index=cells),
    var=pd.DataFrame({
        "gene_biotype": (
            ["protein_coding"] * 2 + ["tRNA"] * 5
        ),
        "is_tRNA": [False] * 2 + [True] * 5,
    }, index=[
        "geneA", "geneB",
        "trna1", "trna2", "trna3", "trna4", "trna5",
    ]),
)
rna.layers["counts"] = csr_matrix(counts)
rna.obsm["X_umap"] = rng.normal(size=(n, 2))

end_counts = np.array(
    [[10, 2]] * 20 + [[2, 10]] * 20
)
end5 = ad.AnnData(
    X=csr_matrix(end_counts),
    obs=pd.DataFrame({
        "present_in_modality": True,
    }, index=qualified),
    var=pd.DataFrame(index=["locusA", "locusB"]),
)

with tempfile.TemporaryDirectory() as directory:
    tmp = Path(directory)
    rna_file = tmp / "cluster.h5ad"
    mu_file = tmp / "sample.h5mu"

    rna.write_h5ad(rna_file)
    mu.MuData({"end5": end5}).write_h5mu(mu_file)

    run(
        "embed-trna.py",
        "--clustered", rna_file,
        "--output", tmp / "trna",
    )
    assert (tmp / "trna.h5ad").exists()
    assert (tmp / "trna.umap.png").exists()
    embedded = ad.read_h5ad(tmp / "trna.h5ad")
    assert embedded.n_obs == 40
    assert embedded.n_vars == 5
    assert "X_umap" in embedded.obsm
    print("PASS: independent tRNA embedding")

    run(
        "plot-end5-usage.py",
        "--clustered", rna_file,
        "--mudata", mu_file,
        "--output", tmp / "end5",
    )
    usage = pd.read_csv(
        tmp / "end5.cluster-usage.csv", index_col=0
    )
    assert np.isclose(usage.loc["locusA", "0"], 10 / 12)
    assert np.isclose(usage.loc["locusA", "1"], 2 / 12)
    assert (tmp / "end5.heatmap.png").exists()
    print("PASS: independent 5'-end usage")

    run(
        "plot-gene-correlations.py",
        "--clustered", rna_file,
        "--output", tmp / "correlations",
    )
    corr = pd.read_csv(
        tmp / "correlations.matrix.csv", index_col=0
    )
    assert np.isclose(corr.loc["geneA", "geneB"], 1)
    assert (tmp / "correlations.heatmap.png").exists()
    print("PASS: independent gene correlations")
