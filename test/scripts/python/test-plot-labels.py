#!/usr/bin/env python
"""Verify readable labels preserve stable identifiers."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "bin/viz"))

from viz_utils import gene_label, end_locus_label

gene_id = "GCF_000005845.2:ASM584v2:NC_000913.3:b0014"

assert gene_label(gene_id, "dnaK", "b0014") == "dnaK (b0014)"
assert gene_label(gene_id, "", "b0014") == "b0014"
assert gene_label(gene_id, "b0014", "b0014") == "b0014"

assert end_locus_label(
    "end5_001",
    "ID=gene-b0014;gene=dnaK;locus_tag=b0014",
) == "dnaK (end5_001)"

assert end_locus_label("end5_002", "") == "end5_002"

# Verify numeric chromosome-position sorting.
loci = pd.DataFrame({
    "contig": ["chr1", "chr1", "chr2", "chr1"],
    "start1": [1000, 100, 50, 200],
    "end1": [1010, 110, 60, 210],
    "_index": [0, 1, 2, 3],
})

ordered = loci.sort_values(
    ["contig", "start1", "end1"],
    kind="stable",
)

assert ordered["_index"].tolist() == [1, 3, 0, 2]

# Relabelling should preserve matrix values and order.
corr = pd.DataFrame(
    [[1.0, 0.4], [0.4, 1.0]],
    index=["geneA", "geneB"],
    columns=["geneA", "geneB"],
)
display = corr.copy()
display.index = ["dnaK (geneA)", "gyrA (geneB)"]
display.columns = display.index

np.testing.assert_array_equal(
    corr.to_numpy(), display.to_numpy()
)

print("PASS: gene labels and genomic end-locus ordering")