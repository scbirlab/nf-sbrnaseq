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


SCRIPT = find_repo_root(Path(__file__), "bin/anndata/build-mudata.py")

def make_object(barcodes, counts, modality=None):
    obj = ad.AnnData(
        X=csr_matrix(np.asarray(counts, dtype=np.int64)),
        obs=pd.DataFrame(index=barcodes),
    )
    if modality:
        obj.uns["__modality__"] = modality
    return obj

with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)

    # RNA has A, B, C; only A and C pass QC.
    rna = make_object(
        ["A", "B", "C"],
        [[1, 2], [3, 4], [5, 6]],
    )
    qc = rna[["A", "C"], :].copy()

    # End loci have A and D; modification sites B and D.
    end5 = make_object(
        ["D", "A"], [[7], [8]], "end5"
    )
    mod = make_object(
        ["B", "D"], [[9], [10]], "mod"
    )

    inputs = {
        "rna": rna,
        "qc": qc,
        "end5": end5,
        "mod": mod,
    }
    for name, obj in inputs.items():
        obj.write_h5ad(tmp / f"{name}.h5ad")

    outfile = tmp / "sample.h5mu"

    subprocess.run([
        sys.executable, str(SCRIPT),
        "--rna", str(tmp / "rna.h5ad"),
        "--rna-qc", str(tmp / "qc.h5ad"),
        "--modalities",
        str(tmp / "end5.h5ad"),
        str(tmp / "mod.h5ad"),
        "--sample-id", "sample1",
        "--output", str(outfile),
    ], check=True)

    result = mu.read_h5mu(outfile)

    expected = ["sample1:A", "sample1:B",
                "sample1:C", "sample1:D"]

    assert list(result.obs_names) == expected

    assert result.obs["rna_qc_pass"].tolist() == [
        True, False, True, False
    ]
    assert result.obs["has_rna"].tolist() == [
        True, True, True, False
    ]
    assert result.obs["has_end5"].tolist() == [
        True, False, False, True
    ]
    assert result.obs["has_mod"].tolist() == [
        False, True, False, True
    ]

    for key, original in (
        ("rna", rna), ("end5", end5), ("mod", mod)
    ):
        actual = result.mod[key]
        assert list(actual.obs_names) == expected
        assert actual.X.sum() == original.X.sum()

    assert result.mod["rna"].X[3, :].sum() == 0
    assert not result.mod["rna"].obs.loc[
        "sample1:D", "present_in_modality"
    ]

    print("PASS: barcode alignment")
    print("PASS: RNA-QC membership")
    print("PASS: missing-modality flags")
    print("PASS: sparse count conservation")
