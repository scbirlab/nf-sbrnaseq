
#!/usr/bin/env python3
"""Assemble RNA, 5'-end and modification AnnData."""

import argparse

import anndata as ad
import mudata as mu
import numpy as np
from scipy.sparse import csr_matrix, issparse


def qualify(mdata, sample_id):
    """Validate barcodes and construct sample-qualified IDs."""
    mdata = mdata.copy()
    barcodes = mdata.obs_names.astype(str)

    if not sample_id or ":" in sample_id:
        raise ValueError("Invalid sample identifier")

    if (
        mdata.obs_names.has_duplicates
        or any(not x or x.lower() in ("nan", "none") for x in barcodes)
    ):
        raise ValueError("Missing or duplicate cell barcodes")

    if "sample_id" in mdata.obs:
        existing = mdata.obs["sample_id"].dropna().astype(str)
        if not existing.eq(sample_id).all():
            raise ValueError("Conflicting sample identifiers")

    if "cell_barcode" in mdata.obs:
        existing = mdata.obs["cell_barcode"].astype(str)
        if not np.array_equal(existing, barcodes):
            raise ValueError("Conflicting cell barcode annotations")

    mdata.obs["cell_barcode"] = list(barcodes)
    mdata.obs["sample_id"] = sample_id
    mdata.obs_names = [f"{sample_id}:{x}" for x in barcodes]
    return mdata



def reindex_sparse(mdata, names, sample_id):
    """Expand rows without densifying X or its layers."""
    mdata = mdata.copy()
    old = mdata.obs_names.get_indexer(names)

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

    X = expand(mdata.X)
    layers = {
        key: expand(value)
        for key, value in mdata.layers.items()
    }

    
    import pandas as pd

    obs = mdata.obs.reindex(names).copy()

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
        var=mdata.var.copy(),
        uns=mdata.uns.copy(),
    )
    for key, value in layers.items():
        output.layers[key] = value

    return output



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rna", required=True)
    parser.add_argument("--rna-qc", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--modalities", nargs="*", default=[])
    args = parser.parse_args()

    rna = qualify(ad.read_h5ad(args.rna), args.sample_id)
    rna_qc = qualify(
        ad.read_h5ad(args.rna_qc), args.sample_id
    )

    if not rna_qc.obs_names.isin(rna.obs_names).all():
        raise ValueError(
            "RNA-QC barcodes are absent from unfiltered RNA"
        )

    objects = {"rna": rna}

    for filename in args.modalities:
        adata = ad.read_h5ad(filename)
        key = str(adata.uns["__modality__"])

        if key in objects:
            raise ValueError(f"Duplicate modality: {key}")

        objects[key] = qualify(adata, args.sample_id)

    # Keep all observed barcodes. Never interpret a
    # structurally absent modality as measured zero.
    names = sorted(set().union(*(
        set(_mdata.obs_names) for _mdata in objects.values()
    )))

    objects = {
        key: reindex_sparse(_mdata, names, args.sample_id)
        for key, _mdata in objects.items()
    }

    # Validate sparse reindexing and record modality presence.
    for key, _mdata in objects.items():
        if list(_mdata.obs_names) != names:
            raise ValueError(f"Misaligned modality: {key}")

        if "present_in_modality" not in _mdata.obs:
            raise ValueError(f"Missing presence mask: {key}")

        present = _mdata.obs["present_in_modality"].to_numpy(dtype=bool)

        absent = _mdata[~present, :]
        if absent.X.nnz if issparse(absent.X) else np.any(absent.X):
            raise ValueError(
                f"Nonzero counts in absent modality: {key}"
            )

    mdata = mu.MuData(objects)
    mdata.uns["observation_join"] = "union"
    mdata.uns["sample_id"] = args.sample_id
    mdata.update()

    if list(mdata.obs_names) != names:
        raise ValueError("MuData observation alignment failed")

    qc_names = set(rna_qc.obs_names)

    mdata.obs["rna_qc_pass"] = [
        name in qc_names for name in names
    ]

    for key, _mdata in objects.items():
        mdata.obs[f"has_{key}"] = (
            _mdata.obs["present_in_modality"]
            .to_numpy(dtype=bool)
        )

    mdata.uns["__modalities__"] = np.asarray(list(objects), dtype=str)
    if (mdata.obs["rna_qc_pass"] > mdata.obs["has_rna"]).any():
        raise ValueError("Invalid RNA-QC membership")

    mdata.write_h5mu(args.output)
    print(mdata)


if __name__ == "__main__":
    main()
