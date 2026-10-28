
"""Reusable clustered heatmaps."""

import matplotlib.pyplot as plt
import seaborn as sns


def heatmap(
    df,
    color_label,
    filename,
    output_format="both",
    row_idx=None,
    col_idx=None,
    *,
    row_colors=None,
    col_colors=None,
    row_linkage=None,
    col_linkage=None,
    cmap="cividis",
    vmin=0.,
    vmax=1.,
    center=None,
    figsize=None,
    mask_missing=True,
):
    """Plot a clustered heatmap with optional external annotations.

    Defaults preserve the original viz_utils heatmap behaviour.
    Explicit linkage matrices permit identical clustering across
    multiple heatmaps of the same underlying data.
    """
    if row_idx is not None:
        df = df.iloc[row_idx]
        if row_colors is not None:
            row_colors = row_colors.reindex(df.index)

    if col_idx is not None:
        df = df.iloc[:, col_idx]
        if col_colors is not None:
            col_colors = col_colors.reindex(df.columns)

    missing = df.isna()
    values = df.fillna(0.)

    row_cluster = (
        row_linkage is not None
        or (
            row_idx is None
            and values.shape[0] > 1
            and values.drop_duplicates().shape[0] > 1
        )
    )
    col_cluster = (
        col_linkage is not None
        or (
            col_idx is None
            and values.shape[1] > 1
            and values.T.drop_duplicates().shape[0] > 1
        )
    )

    if figsize is None:
        figsize = (
            max(5., min(24., values.shape[1] * .2)),
            max(5., min(24., values.shape[0] * .2)),
        )

    hm = sns.clustermap(
        values,
        metric="euclidean",
        row_cluster=row_cluster,
        col_cluster=col_cluster,
        row_linkage=row_linkage,
        col_linkage=col_linkage,
        row_colors=row_colors,
        col_colors=col_colors,
        cbar_kws={"label": color_label},
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        center=center,
        xticklabels=True,
        yticklabels=True,
        mask=missing if mask_missing else None,
        figsize=figsize,
    )

    for _ax in ("x", "y"):
        hm.ax_heatmap.tick_params(
            axis=_ax, 
            labelsize=4.,
        )

    formats = (
        ("png", "pdf")
        if output_format == "both"
        else (output_format,)
    )
    if any(fmt not in ("png", "pdf") for fmt in formats):
        plt.close(hm.fig)
        raise ValueError(
            "output_format must be 'png', 'pdf' or 'both'"
        )

    for fmt in formats:
        hm.savefig(
            f"{filename}.{fmt}",
            dpi=600,
            bbox_inches="tight",
        )

    return hm
