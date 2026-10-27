
import subprocess
import sys
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from test_utils.files import find_repo_root

SCRIPT = find_repo_root(Path(__file__), "bin/anndata-utils.py")

counts = np.array([
    [40, 30, 0],   # 70 UMIs: retain
    [100, 0, 0],   # 100 UMIs: retain
    [0, 0, 5],     # 5 UMIs: exclude
    [30, 35, 0],   # 65 UMIs: retain
], dtype=np.int64)

with tempfile.TemporaryDirectory() as directory:
    directory = Path(directory)
    infile = directory / "input.h5ad"
    prefix = directory / "test"

    adata = ad.AnnData(
        X=counts.copy(),
        obs=pd.DataFrame(
            index=[f"cell_{i}" for i in range(4)]
        ),
        var=pd.DataFrame(
            {"genome_accession": ["test"] * 3},
            index=["gene_a", "gene_b", "gene_c"],
        ),
    )
    adata.layers["counts"] = counts.copy()
    adata.write_h5ad(infile)

    subprocess.run([
        sys.executable,
        str(SCRIPT),
        "filter",
        str(infile),
        "--min-counts-per-cell", "60",
        "--min-cells-per-gene", "2",
        "-o", str(prefix),
    ], check=True, cwd=directory)

    result = ad.read_h5ad(
        f"{prefix}.filtered.h5ad"
    )

    assert result.n_obs == 3
    assert result.n_vars == 2

    assert set(result.obs_names) == {
        "cell_0", "cell_1", "cell_3"
    }
    assert set(result.var_names) == {
        "gene_a", "gene_b"
    }

    np.testing.assert_array_equal(
        result.X,
        result.layers["counts"],
    )

    assert result.uns["filtering"]["min_counts_per_cell"] == 60
    assert result.uns["filtering"]["min_cells_per_gene"] == 2

    retention = pd.read_csv(
        f"{prefix}.cell-retention.csv"
    )
    assert retention.loc[
        retention.min_counts_per_cell == 60,
        "n_cells"
    ].iloc[0] == 3

    print("PASS: explicit thresholds")
    print("PASS: cell and gene retention")
    print("PASS: raw counts preserved")
    print("PASS: retention curve")
