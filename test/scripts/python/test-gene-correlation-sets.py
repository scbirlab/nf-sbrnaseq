
#!/usr/bin/env python
"""Synthetic regression test for three annotation heatmaps."""

import subprocess
import sys
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "bin/viz/plot-gene-correlations.py"

with tempfile.TemporaryDirectory() as temp:
    directory = Path(temp)
    rng = np.random.default_rng(42)

    counts = rng.poisson(5, size=(40, 12)).astype(float)
    genes = [
        f"acc:asm:chr1:g{i}" for i in range(12)
    ]

    adata = ad.AnnData(
        X=csr_matrix(np.log1p(counts)),
        obs=pd.DataFrame(
            index=[f"cell{i}" for i in range(40)]
        ),
        var=pd.DataFrame(index=genes),
    )
    adata.layers["counts"] = csr_matrix(counts)

    h5ad = directory / "rna.h5ad"
    adata.write_h5ad(h5ad)

    annotations = pd.DataFrame({
        "genome_accession": ["acc"] * 12,
        "assembly": ["asm"] * 12,
        "chr": ["chr1"] * 12,
        "gene_id": [f"g{i}" for i in range(12)],
        "COG_category": [
            "JK", "J", "E", "E", "F", "F",
            "L", "L", "S", "S", "", "",
        ],
        "EC": [
            "1.1.1.1,2.1.1.1", "1.1.1.1",
            "2.1.1.1", "2.1.1.1",
            "3.1.1.1", "3.1.1.1",
            "", "", "", "", "", "",
        ],
        "KEGG_Pathway": [
            "ko00010,ko00020", "ko00010",
            "ko00020", "ko00020",
            "ko00030", "ko00030",
            "", "", "", "", "", "",
        ],
    })
    annotations["KEGG_Module"] = ""
    annotations["BRITE"] = ""

    table = directory / "functional.tsv"
    annotations.to_csv(table, sep="\t", index=False)

    prefix = directory / "correlations"

    subprocess.run([
        sys.executable, str(SCRIPT),
        "--clustered", str(h5ad),
        "--annotations", str(table),
        "--output", str(prefix),
    ], check=True)

    for kind in ("cog", "ec", "kegg"):
        assert (
            directory / f"correlations.{kind}.heatmap.png"
        ).exists()
        assert (
            directory / f"correlations.{kind}.colors.csv"
        ).exists()

    metadata = pd.read_csv(
        directory / "correlations.annotations.csv",
        index_col=0,
        keep_default_na=False,
    )

    assert metadata.loc[genes[0], "cog_memberships"] == "J,K"
    assert metadata.loc[genes[0], "cog_group"] == "J"
    for kind in ("cog", "ec", "kegg"):
        assert (
            metadata.loc[genes[11], f"{kind}_heatmap_color"]
            == "lightgrey"
        )
        assert (
            metadata.loc[genes[11], f"{kind}_umap_color"]
            == "lightgrey"
        )

    matrix = pd.read_csv(
        directory / "correlations.matrix.csv",
        index_col=0,
    )
    assert matrix.shape == (12, 12)
    np.testing.assert_allclose(
        np.diag(matrix), np.ones(12)
    )
    np.testing.assert_allclose(
        matrix.to_numpy(), matrix.to_numpy().T
    )

print("PASS: COG, EC and KEGG correlation heatmaps")
