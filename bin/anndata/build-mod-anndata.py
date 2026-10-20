
#!/usr/bin/env python
"""Build sparse cell × modification-site AnnData."""

import argparse
import gzip
import pandas as pd
import numpy as np
import anndata as ad
from scipy.sparse import coo_matrix


COUNT_COLUMNS = (
    "count_A", "count_C", "count_G", "count_T",
    "ref_count", "alt_count", "other_count",
    "callable_depth",
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--modality", default=None, type=str)
    args = parser.parse_args()

    df = pd.read_csv(
        args.counts,
        sep="\t",
        compression="infer",
        low_memory=False,
    )

    required = {
        "sample_id", "cell_barcode", "site_id",
        *COUNT_COLUMNS,
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")

    df = df.loc[df.sample_id == args.sample_id].copy()

    if df.duplicated(["cell_barcode", "site_id"]).any():
        raise ValueError("Duplicate cell-site count rows")

    for column in COUNT_COLUMNS:
        values = pd.to_numeric(df[column], errors="raise")
        if values.isna().any() or (values < 0).any():
            raise ValueError(f"Invalid {column}")
        if (values != np.floor(values)).any():
            raise ValueError(f"Noninteger {column}")
        df[column] = values.astype(np.int64)

    if not (
        df["callable_depth"] == df[[f"count_{b}" for b in "ACGT"]].sum(axis=1)
    ).all():
        raise ValueError("Base-count conservation failed")

    if not (
        (df["callable_depth"] == df["ref_count"]) + df["alt_count"]
        + df["other_count"]
    ).all():
        raise ValueError("REF/ALT count conservation failed")

    cells = pd.Index(
        sorted(df.cell_barcode.astype(str).unique()),
        name="cell_barcode",
    )
    sites = pd.Index(
        sorted(df.site_id.astype(str).unique()),
        name="site_id",
    )

    cell_codes = cells.get_indexer(
        df.cell_barcode.astype(str)
    )
    site_codes = sites.get_indexer(
        df.site_id.astype(str)
    )
    shape = (len(cells), len(sites))

    def sparse_column(column):
        return coo_matrix(
            (
                df[column].to_numpy(dtype=np.int64),
                (cell_codes, site_codes),
            ),
            shape=shape,
            dtype=np.int64,
        ).tocsr()

    ref = sparse_column("ref_count")
    alt = sparse_column("alt_count")
    other = sparse_column("other_count")
    depth = sparse_column("callable_depth")

    # Metadata must be consistent for each resolved site.
    measurement_columns = set(COUNT_COLUMNS) | {
        "sample_id", "cell_barcode", "alt_fraction",
    }
    metadata_columns = [
        c for c in df.columns
        if c not in measurement_columns
    ]
    metadata = df[metadata_columns].drop_duplicates()

    if metadata["site_id"].duplicated().any():
        raise ValueError(
            "Inconsistent metadata for a modification site"
        )

    metadata = (
        metadata.set_index("site_id")
        .reindex(sites)
    )

    adata = ad.AnnData(
        X=alt,
        obs=pd.DataFrame(index=cells),
        var=metadata,
    )
    adata.layers["ref_count"] = ref
    adata.layers["alt_count"] = alt.copy()
    adata.layers["other_count"] = other
    adata.layers["callable_depth"] = depth
    adata.layers["ref_alt_depth"] = ref + alt

    for base in "ACGT":
        adata.layers[f"count_{base}"] = sparse_column(
            f"count_{base}"
        )

    adata.obs["sample_id"] = args.sample_id
    adata.uns["alt_fraction_denominator"] = "ref_alt_depth"

    assert int(adata.X.sum()) == int(df.alt_count.sum())
    assert int(adata.layers["callable_depth"].sum()) == int(df.callable_depth.sum())

    adata.uns["__modality__"] = args.modality or "RNA modification"
    adata.write_h5ad(args.output, compression="gzip")
    print(adata)


if __name__ == "__main__":
    main()
