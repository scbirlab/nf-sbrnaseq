
#!/usr/bin/env python
"""Build sparse cell-by-5'-locus AnnData."""

import argparse

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts", required=True)
    parser.add_argument("--loci", required=True)
    parser.add_argument("--qc", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--modality", default=None, type=str)
    args = parser.parse_args()

    loci = pd.read_csv(args.loci, sep="\t")
    counts = pd.read_csv(args.counts, sep="\t")
    qc = pd.read_csv(args.qc, sep="\t")

    required_loci = {
        "end5_locus_id", "contig", "start1", "end1", "strand"
    }
    missing = required_loci - set(loci.columns)
    if missing:
        raise ValueError(f"Missing locus annotations: {sorted(missing)}")

    if loci["end5_locus_id"].isna().any():
        raise ValueError("Missing locus identifiers")

    loci = loci.set_index("end5_locus_id", verify_integrity=True)

    if (loci.index.astype(str).str.len() == 0).any():
        raise ValueError("Empty locus identifiers")

    if not loci["strand"].isin(["+", "-"]).all():
        raise ValueError("Invalid locus strand")

    start = pd.to_numeric(loci["start1"], errors="raise")
    end = pd.to_numeric(loci["end1"], errors="raise")
    if (start < 1).any() or (end < start).any():
        raise ValueError("Invalid locus coordinates")

    counts = counts.loc[
        counts["sample_id"] == args.sample_id
    ].copy()
    qc = qc.loc[
        qc["sample_id"] == args.sample_id
    ].copy()

    if counts.duplicated(
        ["cell_barcode", "end5_locus_id"]
    ).any():
        raise ValueError("Duplicate cell-locus counts")

    if qc["cell_barcode"].duplicated().any():
        raise ValueError("Duplicate cell QC identifiers")

    unknown = set(counts["end5_locus_id"]) - set(loci.index)
    if unknown:
        raise ValueError(f"Unknown end5 loci: {unknown}")

    if not counts["cell_barcode"].isin(qc["cell_barcode"]).all():
        raise ValueError("Counts contain cells absent from QC")

    if counts["molecule_count"].isna().any():
        raise ValueError("Missing molecule counts")

    values = counts["molecule_count"].to_numpy()
    if np.any(values < 0) or np.any(values != np.floor(values)):
        raise ValueError("Molecule counts must be nonnegative integers")

    cells = pd.Index(
        sorted(qc["cell_barcode"].astype(str)),
        name="cell_barcode",
    )
    features = pd.Index(
        loci.index.astype(str),
        name="end5_locus_id",
    )

    row_lookup = pd.Series(
        np.arange(len(cells)), index=cells
    )
    col_lookup = pd.Series(
        np.arange(len(features)), index=features
    )

    row = row_lookup.loc[
        counts["cell_barcode"].astype(str)
    ].to_numpy()
    col = col_lookup.loc[
        counts["end5_locus_id"].astype(str)
    ].to_numpy()

    matrix = coo_matrix(
        (values.astype(np.int64), (row, col)),
        shape=(len(cells), len(features)),
        dtype=np.int64,
    ).tocsr()

        # Independently validate each cell against assignment QC.
    qc_by_cell = qc.set_index("cell_barcode")

    required_qc = {
        "total_ends_count",
        "assigned_ends_count",
        "unassigned_ends_count",
        "assigned_ends_fraction",
    }
    missing = required_qc - set(qc_by_cell.columns)
    if missing:
        raise ValueError(f"Missing assignment QC: {sorted(missing)}")

    for column in (
        "total_ends_count",
        "assigned_ends_count",
        "unassigned_ends_count",
    ):
        qc_by_cell[column] = pd.to_numeric(
            qc_by_cell[column], errors="raise"
        )
        if (
            qc_by_cell[column].isna().any()
            or (qc_by_cell[column] < 0).any()
            or (qc_by_cell[column] % 1 != 0).any()
        ):
            raise ValueError(f"Invalid {column}")

    total = qc_by_cell["total_ends_count"]
    assigned = qc_by_cell["assigned_ends_count"]
    unassigned = qc_by_cell["unassigned_ends_count"]

    if not (total == assigned + unassigned).all():
        raise ValueError("Assignment QC conservation failed")

    observed = pd.Series(
        np.asarray(matrix.sum(axis=1)).ravel(),
        index=cells,
    )
    expected = assigned.reindex(cells)

    if not np.array_equal(observed.to_numpy(), expected.to_numpy()):
        raise ValueError("Per-cell assigned-end counts disagree with QC")

    fractions = pd.to_numeric(
        qc_by_cell["assigned_ends_fraction"],
        errors="raise",
    )
    expected_fractions = assigned.div(total.where(total > 0))

    if not np.allclose(
        fractions.to_numpy(),
        expected_fractions.to_numpy(),
        rtol=1e-6,
        atol=1e-8,
        equal_nan=True,
    ):
        raise ValueError("Assignment fractions disagree with QC")

    obs = qc.set_index("cell_barcode").reindex(cells)
    for column in (
        "total_ends_count",
        "assigned_ends_count",
        "unassigned_ends_count",
        "assigned_ends_fraction",
    ):
        if column in obs:
            obs[column] = pd.to_numeric(
                obs[column], 
                errors="raise",
            )

    obs["sample_id"] = args.sample_id

    adata = ad.AnnData(
        X=matrix,
        obs=obs,
        var=loci.reindex(features),
    )
    adata.obs_names = cells.astype(str)
    adata.var_names = features.astype(str)

    if int(adata.X.sum()) != int(values.sum()):
        raise ValueError("Count conservation failed")

    adata.uns["__modality__"] = args.modality or "end5"
    adata.write_h5ad(args.output, compression="gzip")

    print(adata)


if __name__ == "__main__":
    main()
