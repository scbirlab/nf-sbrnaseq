
import numpy as np
import scanpy as sc
import anndata as ad

# Two cells with identical transcript proportions
# but tenfold different sequencing depths.
counts = np.array([
    [10, 20, 30],
    [100, 200, 300],
], dtype=np.int64)

adata = ad.AnnData(counts.copy())
adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Raw counts must remain unchanged.
np.testing.assert_array_equal(
    adata.layers["counts"], counts
)

# Both cells should now have identical expression.
np.testing.assert_allclose(
    adata.X[0],
    adata.X[1],
    rtol=1e-6,
    atol=1e-6,
)

# Before log1p, normalised counts must sum to 10,000.
np.testing.assert_allclose(
    np.expm1(adata.X).sum(axis=1),
    [10000, 10000],
    rtol=1e-6,
)

print("PASS: normalisation preserves raw counts")
print("PASS: identical proportions give identical expression")
print("PASS: normalised totals equal 10,000")
