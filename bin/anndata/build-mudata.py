
#!/usr/bin/env python3
"""Assemble RNA, 5'-end and modification AnnData."""

import argparse

import anndata as ad
import mudata as mu
import numpy as np
from scipy.sparse import csr_matrix, issparse


def qualify(obj, sample_id):
    """Create unique sample-qualified observation IDs."""
    obj = obj.copy()
    names = [f"{sample_id}:{x}" for x in obj.obs_names]

    if len(names) != len(set(names)):
        raise ValueError("Duplicate cell identifiers")

    obj.obs["cell_barcode"] = list(obj.obs_names)
    obj.obs["sample_id"] = sample_id
    obj.obs_names = names
    return obj


def reindex_sparse(obj, names, sample_id):
    """Expand rows without densifying X or its layers."""
    obj = obj.copy()
    old = obj.obs_names.get_indexer(names)

    present = old >= 0
    valid_rows = np.flatnonzero(present)
    source_rows = old[present]

    def expand(matrix):
        if not issparse(matrix):
            matrix = csr_matrix(matrix)
        matrix = matrix.tocsr()
        result = csr_matrix(
            (
                matrix[source_rows].data,
                matrix[source_rows].indices,
                matrix[source_rows].indptr,
            ),
            shape=(len(source_rows), matrix.shape[1]),
        )
        # Allocate positions in the common row index.
        from scipy.sparse import coo_matrix
        selected = result.tocoo()
        return coo_matrix(
            (
                selected.data,
                (
                    valid_rows[selected.row],
                    selected.col,
                ),
            ),
            shape=(len(names), matrix.shape[1]),
        ).tocsr()

    X = expand(obj.X)
    layers = {
        key: expand(value)
        for key, value in obj.layers.items()
    }

    
    import pandas as pd

    obs = obj.obs.reindex(names).copy()

    # Restore identifiers for all cells, including cells
    # structurally absent from this modality.
    sample_id = str(sample_id)
    obs["sample_id"] = sample_id
    obs["cell_barcode"] = [
        name[len(sample_id) + 1:]
        for name in names
    ]

    # Preserve missing modality measurements as NaN,
    # but normalise dtypes for HDF5 serialisation.
    for column in obs.columns:
        if column in ("sample_id", "cell_barcode"):
            continue

        series = obs[column]

        if pd.api.types.is_categorical_dtype(series):
            continue

        if pd.api.types.is_numeric_dtype(series):
            obs[column] = pd.to_numeric(
                series, errors="raise"
            )
            continue

        # Object columns may contain strings, numeric
        # values and NaNs after union reindexing.
        nonmissing = series.dropna()

        if len(nonmissing) == 0:
            obs[column] = pd.Series(
                float("nan"), index=obs.index
            )
        elif nonmissing.map(
            lambda x: isinstance(
                x, (int, float, np.integer, np.floating)
            )
        ).all():
            obs[column] = pd.to_numeric(
                series, errors="raise"
            )
        else:
            obs[column] = (
                series.astype("string")
                .fillna("")
                .astype(str)
            )

    obs["present_in_modality"] = present



    output = ad.AnnData(
        X=X,
        obs=obs,
        var=obj.var.copy(),
        uns=obj.uns.copy(),
    )
    for key, value in layers.items():
        output.layers[key] = value

    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rna", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--modalities", nargs='*')
    args = parser.parse_args()

    objects = {
        "rna": qualify(
            ad.read_h5ad(args.rna), args.sample_id
        )
    }

    for f in args.modalities:
        adata = ad.read_h5ad(f)
        key = adata.uns["__modality__"]
        objects[key] = qualify(adata, args.sample_id)

    names = sorted(set().union(*(
        set(obj.obs_names)
        for obj in objects.values()
    )))

    objects = {
        key: reindex_sparse(obj, names, args.sample_id)
        for key, obj in objects.items()
    }

    mdata = mu.MuData(objects)
    mdata.uns["observation_join"] = "union"
    mdata.uns["sample_id"] = args.sample_id
    mdata.update()
    mdata.write_h5mu(args.output)
    print(mdata)


if __name__ == "__main__":
    main()
