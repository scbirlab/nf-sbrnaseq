
#!/usr/bin/env python
"""Regression test for functional heatmaps and gene UMAPs."""

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

    genes = [f"acc:asm:chr1:g{i}" for i in range(20)]
    counts = rng.poisson(5, size=(40, 20)).astype(float)

    rna = ad.AnnData(
        X=csr_matrix(np.log1p(counts)),
        obs=pd.DataFrame(
            index=[f"cell{i}" for i in range(40)]
        ),
        var=pd.DataFrame(index=genes),
    )
    rna.layers["counts"] = csr_matrix(counts)

    h5ad = directory / "rna.h5ad"
    rna.write_h5ad(h5ad)

    annotations = pd.DataFrame({
        "genome_accession": ["acc"] * 20,
        "assembly": ["asm"] * 20,
        "chr": ["chr1"] * 20,
        "gene_id": [f"g{i}" for i in range(20)],
        "COG_category": [
            "JK", "J", "E", "E", "F", "F",
            "L", "L", "S", "S", "O", "O",
            "G", "G", "C", "C", "M", "M", "", "",
        ],
        "EC": [
            f"{i % 10 + 1}.1.1.1" for i in range(18)
        ] + ["", ""],
        "KEGG_Pathway": [
            f"ko{i % 10:05d}" for i in range(18)
        ] + ["", ""],
    })

    table = directory / "functional.tsv"
    annotations.to_csv(table, sep="\t", index=False)

    prefix = directory / "correlations"

    subprocess.run([
        sys.executable, str(SCRIPT),
        "--clustered", str(h5ad),
        "--annotations", str(table),
        "--output", str(prefix),
        "--top-groups", "8",
    ], check=True)

    for kind in ("cog", "ec", "kegg"):
        assert (
            directory / f"correlations.{kind}.heatmap.png"
        ).is_file()
        assert (
            directory / f"correlations.{kind}.gene-umap.png"
        ).is_file()
        assert (
            directory / f"correlations.{kind}.colors.csv"
        ).is_file()

        colors = pd.read_csv(
            directory / f"correlations.{kind}.colors.csv"
        )
        assert colors["highlighted_on_umap"].sum() <= 8

    metadata = pd.read_csv(
        directory / "correlations.annotations.csv",
        index_col=0,
        keep_default_na=False,
    )
    coords = pd.read_csv(
        directory / "correlations.gene-umap.csv",
        index_col=0,
    )

    assert len(metadata) == 20
    assert metadata.loc[genes[0], "cog_memberships"] == "J,K"
    assert metadata.loc[genes[0], "cog_group"] == "J"

    for kind in ("cog", "ec", "kegg"):
        assert (
            metadata.loc[genes[19], f"{kind}_umap_color"]
            == "lightgrey"
        )

    np.testing.assert_allclose(
        metadata.loc[coords.index, ["UMAP_1", "UMAP_2"]]
        .to_numpy(dtype=float),
        coords.to_numpy(dtype=float),
    )
    assert np.isfinite(coords.to_numpy()).all()

    matrix = pd.read_csv(
        directory / "correlations.matrix.csv",
        index_col=0,
    )
    assert matrix.shape == (20, 20)
    np.testing.assert_allclose(
        matrix.to_numpy(),
        matrix.to_numpy().T,
        atol=1e-12,
    )

print("PASS: three heatmaps, three gene UMAPs, shared coordinates")
