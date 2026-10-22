import seaborn as sns


def heatmap(
    df, 
    color_label, 
    filename,
    output_format="png",
    row_idx=None, 
    col_idx=None
):
    if row_idx is not None:
        df = df.iloc[row_idx]
    if col_idx is not None:
        df = df.iloc[:,col_idx]
    missing = df.isna()
    values = df.fillna(0.)
    # Only cluster an axis containing distinguishable observations.
    row_cluster = (
        row_idx is None
        and values.shape[0] > 1
        and values.drop_duplicates().shape[0] > 1
    )
    col_cluster = (
        col_idx is None
        and values.shape[1] > 1
        and values.T.drop_duplicates().shape[0] > 1
    )
    hm = sns.clustermap(
        values,
        metric="euclidean",
        row_cluster=row_cluster,
        col_cluster=col_cluster,
        cbar_kws={"label": color_label},
        cmap="cividis",
        vmin=0., 
        vmax=1.,
        xticklabels=True,
        yticklabels=True,
        mask=missing,
        figsize=(
            max(6, min(14, values.shape[1] * .5)),
            max(6, min(14, values.shape[0] * .25)),
        ),
    )
    hm.ax_heatmap.tick_params(axis="y", labelsize=6.)
    hm.savefig(
        f"{filename}.png", 
        dpi=300, 
        bbox_inches="tight",
    )
    hm.savefig(
        f"{filename}.pdf", 
        dpi=300, 
        bbox_inches="tight",
    )
    return hm