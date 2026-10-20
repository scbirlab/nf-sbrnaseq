
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

    loci = loci.set_index("locus_id", verify_integrity=True)

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

    obs = qc.set_index("cell_barcode").reindex(cells)
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

    adata.uns["__modality__"] = args.modality or "5'-end"
    adata.write_h5ad(args.output, compression="gzip")

    print(adata)


if __name__ == "__main__":
    main()
