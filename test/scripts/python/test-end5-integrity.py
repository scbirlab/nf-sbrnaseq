
#!/usr/bin/env python
"""End-to-end validation of 5'-end assignments and AnnData."""

import gzip
import subprocess
import sys
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from test_utils.files import find_repo_root


ASSIGN = find_repo_root(Path(__file__), "bin/assign-end-loci.py")
BUILD = find_repo_root(Path(__file__), "bin/anndata/build-end5-anndata.py")
MUDATA = find_repo_root(Path(__file__), "bin/anndata/build-mudata.py")


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


with tempfile.TemporaryDirectory() as directory:
    tmp = Path(directory)

    loci = pd.DataFrame([
        {
            "end5_locus_id": "plus",
            "contig": "chr1",
            "strand": "+",
            "start1": 100,
            "end1": 105,
            "support": 3,
        },
        {
            "end5_locus_id": "minus",
            "contig": "chr1",
            "strand": "-",
            "start1": 200,
            "end1": 205,
            "support": 2,
        },
    ])

    # Molecule-level coordinates are 1-based.
    ends = pd.DataFrame([
        ("sample1", "A", "mol1", "chr1", "+", 100),
        ("sample1", "A", "mol2", "chr1", "+", 105),
        ("sample1", "B", "mol3", "chr1", "-", 200),
        ("sample1", "B", "mol4", "chr1", "-", 205),
        ("sample1", "C", "mol5", "chr1", "+", 300),
    ], columns=[
        "sample_id", "cell_barcode", "molecule_id",
        "contig", "strand", "pos1",
    ])

    raw_loci = tmp / "raw-loci.tsv"
    annotation = tmp / "annotated.tsv"
    end_file = tmp / "ends.tsv.gz"
    counts_file = tmp / "counts.tsv"
    qc_file = tmp / "assignment-qc.tsv"
    output = tmp / "end5.h5ad"

    loci.to_csv(raw_loci, sep="\t", index=False)
    loci.to_csv(annotation, sep="\t", index=False)

    with gzip.open(end_file, "wt") as handle:
        ends.to_csv(handle, sep="\t", index=False)

    run(
        ASSIGN,
        "--loci", raw_loci,
        "--ends", end_file,
        "--counts", counts_file,
        "--qc", qc_file,
    )

    counts = pd.read_csv(counts_file, sep="\t")
    qc = pd.read_csv(qc_file, sep="\t")

    assert counts["molecule_count"].sum() == 4
    assert qc["total_ends_count"].sum() == 5
    assert qc["assigned_ends_count"].sum() == 4
    assert qc["unassigned_ends_count"].sum() == 1

    run(
        BUILD,
        "--counts", counts_file,
        "--loci", annotation,
        "--qc", qc_file,
        "--sample-id", "sample1",
        "--output", output,
    )

    obj = ad.read_h5ad(output)

    assert obj.n_obs == 3
    assert obj.n_vars == 2
    assert int(obj.X.sum()) == 4

    assert obj.X[obj.obs_names.get_loc("A"), 
                 obj.var_names.get_loc("plus")] == 2
    assert obj.X[obj.obs_names.get_loc("B"),
                 obj.var_names.get_loc("minus")] == 2
    assert obj.X[obj.obs_names.get_loc("C"), :].sum() == 0

    assert obj.var.loc["plus", "start1"] == 100
    assert obj.var.loc["minus", "strand"] == "-"

    print("PASS: positive/negative strand assignment")
    print("PASS: inclusive locus boundaries")
    print("PASS: assigned/unassigned molecule conservation")
    print("PASS: zero-assignment cell retained")
    print("PASS: locus annotations preserved")

    # Corrupt an individual cell's assigned count.
    invalid_qc = qc.copy()
    invalid_qc.loc[
        invalid_qc.cell_barcode == "A",
        "assigned_ends_count",
    ] += 1
    bad_qc = tmp / "bad-qc.tsv"
    invalid_qc.to_csv(bad_qc, sep="\t", index=False)

    result = run(
        BUILD,
        "--counts", counts_file,
        "--loci", annotation,
        "--qc", bad_qc,
        "--sample-id", "sample1",
        "--output", tmp / "invalid.h5ad",
        success=False,
    )
    assert (
        "conservation failed" in result.stderr
        or "disagree with QC" in result.stderr
    )

    print("PASS: inconsistent assignment QC rejected")

    # Duplicate locus annotations must fail.
    duplicate_loci = pd.concat(
        [loci, loci.iloc[[0]]], ignore_index=True
    )
    bad_loci = tmp / "duplicate-loci.tsv"
    duplicate_loci.to_csv(
        bad_loci, sep="\t", index=False
    )

    result = run(
        BUILD,
        "--counts", counts_file,
        "--loci", bad_loci,
        "--qc", qc_file,
        "--sample-id", "sample1",
        "--output", tmp / "invalid.h5ad",
        success=False,
    )
    assert result.returncode != 0

    print("PASS: duplicate locus annotations rejected")


    # Confirm the 5'-end modality survives MuData assembly.
    import mudata as mu
    from scipy.sparse import csr_matrix

    rna = ad.AnnData(
        X=csr_matrix([[5, 2], [3, 1]]),
        obs=pd.DataFrame(index=["A", "D"]),
        var=pd.DataFrame(index=["gene1", "gene2"]),
    )
    rna.layers["counts"] = rna.X.copy()

    rna_file = tmp / "rna.h5ad"
    rna_qc_file = tmp / "rna-qc.h5ad"
    mudata_file = tmp / "sample.h5mu"

    rna.write_h5ad(rna_file)
    rna[["A"], :].copy().write_h5ad(rna_qc_file)

    run(
        MUDATA,
        "--rna", rna_file,
        "--rna-qc", rna_qc_file,
        "--modalities", output,
        "--sample-id", "sample1",
        "--output", mudata_file,
    )

    mdata = mu.read_h5mu(mudata_file)
    joined = mdata.mod["end5"]

    assert list(mdata.obs_names) == [
        "sample1:A", "sample1:B",
        "sample1:C", "sample1:D",
    ]
    assert int(joined.X.sum()) == 4
    assert joined.var.loc["plus", "start1"] == 100
    assert joined.var.loc["minus", "strand"] == "-"
    assert mdata.obs.loc["sample1:A", "rna_qc_pass"]
    assert not mdata.obs.loc["sample1:D", "has_end5"]
    assert joined["sample1:D", :].X.sum() == 0

    print("PASS: MuData count and annotation preservation")
    print("PASS: structurally absent 5'-end modality")
