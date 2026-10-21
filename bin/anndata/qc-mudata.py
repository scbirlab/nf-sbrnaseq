#!/usr/bin/env python
"""Generate preliminary multimodal QC figures from MuData."""

import argparse
from pathlib import Path

from carabiner.mpl import grid, figsaver
import matplotlib.pyplot as plt
import mudata as mu
import numpy as np
import pandas as pd


def totals(matrix):
    return np.asarray(matrix.sum(axis=1)).ravel()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("--output", default="qc")
    parser.add_argument("--plot-format", default="png")
    args = parser.parse_args()

    _figsaver = figsaver(format=args.plot_format)

    mdata = mu.read_h5mu(args.input)
    prefix = Path(args.output)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    obs = mdata.obs.copy()
    if "rna_qc_pass" not in obs:
        raise ValueError("Missing rna_qc_pass")

    flags = [f"has_{key}" for key in mdata.mod]
    for flag in flags:
        if flag not in obs:
            raise ValueError(f"Missing {flag}")

    # Explicitly retain all union barcodes for QC.
    qc = pd.DataFrame(index=mdata.obs_names)
    qc["rna_qc_pass"] = obs["rna_qc_pass"].astype(bool)

    for flag in flags:
        qc[flag] = obs[flag].astype(bool)

    for key, obj in mdata.mod.items():
        matrix = (
            obj.layers["counts"]
            if key == "rna" and "counts" in obj.layers
            else obj.layers["callable_depth"]
            if key == "mod" and "callable_depth" in obj.layers
            else obj.X
        )
        qc[f"{key}_depth"] = totals(matrix)

    # Modality overlap and cell retention.
    summary = pd.DataFrame([
        {"metric": "union_barcodes", "value": len(qc)},
        {"metric": "rna_qc_pass", "value": int(qc.rna_qc_pass.sum())},
        *[
            {"metric": flag, "value": int(qc[flag].sum())}
            for flag in flags
        ],
    ])

    fig, ax = grid()
    ax.bar(
        "metric", "value",
        data=summary,
    )
    ax.tick_params(axis="x", rotation=40.)
    ax.set(
        ylabel="Number of barcodes",
    )
    _figsaver(
        fig=fig,
        name=f"{prefix}.barcodes",
        df=summary,
    )

    # Correlations between total modality depths.
    depths = [f"{key}_depth" for key in mdata.mod]
    if len(depths) > 1:
        fig, axes = grid(
            ncol=len(depths) - 1,
            squeeze=False,
        )
        for ax, column in zip(fig.axes, depths[1:]):
            ax.scatter(
                "rna_depth", 
                column,
                data=qc,
                s=3., 
                alpha=0.35,
            )
            ax.set(
                xlabel="RNA counts",
                ylabel=column,
                xscale="symlog",
                yscale="symlog",
            )
        _figsaver(
            fig=fig,
            name=f"{prefix}.depths",
            df=qc,
        )

    if "mod" in mdata.mod:
        mod = mdata.mod["mod"]
        depth = mod.layers["callable_depth"]
        alt = mod.layers["alt_count"]
        denom = mod.layers["ref_alt_depth"]

        site_depth = np.asarray(depth.sum(axis=0)).ravel()
        site_alt = np.asarray(alt.sum(axis=0)).ravel()
        site_denom = np.asarray(denom.sum(axis=0)).ravel()

        sites = mod.var.copy()
        sites["callable_depth"] = site_depth
        sites["alt_count"] = site_alt
        sites["ref_alt_depth"] = site_denom
        sites["alt_fraction"] = np.divide(
            site_alt,
            site_denom,
            out=np.full(len(site_alt), np.nan),
            where=site_denom > 0,
        )
        sites.to_csv(f"{prefix}.mod-sites.csv")

        fig, axes = grid(ncol=2)
        axes[0].hist(
            np.log1p(sites["callable_depth"]),
            bins=np.geomspace(
                max(1, np.log1p(sites["callable_depth"]).min()), 
                max(2, np.log1p(sites["callable_depth"]).max()), 
                num=40,
            ),
        )
        axes[0].set(
            xlabel="log1p(callable depth)",
            ylabel="Number of sites",
        )
        axes[1].hist(
            sites["alt_fraction"].dropna(), 
            bins=30, 
            range=(0, 1))
        axes[1].set(
            xlabel="ALT / (REF + ALT)",
            ylabel="Number of sites",
        )
        _figsaver(
            fig=fig,
            name=f"{prefix}.modifications",
            df=sites,
        )

    if "end5" in mdata.mod:
        end5 = mdata.mod["end5"]
        cells = end5.obs

        if "assigned_ends_fraction" in cells:
            fractions = pd.to_numeric(
                cells.loc[
                    cells.present_in_modality.astype(bool),
                    "assigned_ends_fraction"
                ],
                errors="coerce",
            ).dropna()

            fig, ax = grid()
            ax.hist(fractions, bins=30, range=(0, 1))
            ax.set(
                xlabel="Fraction of ends assigned to loci",
                ylabel="Number of barcodes",
            )
            _figsaver(
                fig=fig,
                name=f"{prefix}.end5-assignment",
                df=cells,
            )

        end5.var.to_csv(f"{prefix}.end5-loci.csv")

    print(f"Generated QC report: {prefix}")


if __name__ == "__main__":
    main()
