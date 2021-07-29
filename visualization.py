from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr, spearmanr


def scatter_pred_vs_true(
    x,
    y,
    color: str,
    title: str,
    xlab: str = "Experimental Stability Change",
    ylab: str = "Predicted Stability Change",
    ylim: Tuple[float] = [-22.0, 22.0],
    scatter_size: float = 0.5,
    fontsize: float = 15.0,
    labelpad: float = 10.0,
    title_gap: float = 1.03,
):
    """
    Function to scatter plot true vs predicted ddG and report Pearson's and
    Spearman's correlations.
    """
    fig, ax = plt.subplots()
    ax.scatter(x, y, s=scatter_size, c=color)

    ax.set_ylim(ylim)

    # Make x axis symmetric around 0
    max_abs_lim = max([abs(x) for x in ax.get_xlim()])
    ax.set_xlim([-max_abs_lim, max_abs_lim])

    ax.set_xlabel(xlab, fontsize=fontsize, labelpad=labelpad)
    ax.set_ylabel(ylab, fontsize=fontsize, labelpad=labelpad)

    x_ = [val[0] for val in zip(x, y) if not np.isnan(val[1]) and not np.isnan(val[0])]
    y_ = [val[1] for val in zip(x, y) if not np.isnan(val[1]) and not np.isnan(val[0])]

    title = (
        title
        + f"\nPearson's r {pearsonr(x_, y_)[0]:4.2f}"
        # + f"\nPearson's r {pearsonr(x_, y_)[0]:4.2f}, Spearman's r {spearmanr(x_,y_)[0]:4.2f}"
    )

    ax.set_title(title, fontsize=fontsize, y=title_gap)

    fig.tight_layout()

    return fig, ax


def plot_validation_performance(fig_title, results_dict):
    colors_dict = {
        "dms": "steelblue",
        "protein_g": "firebrick",
        "guerois": "forestgreen",
    }

    fontsize = 15
    fig, ax = plt.subplots()

    for val_key in results_dict.keys():
        pearson_list = results_dict[val_key]
        ax.plot(
            np.arange(len(pearson_list)),
            pearson_list,
            label=val_key,
            color=colors_dict[val_key],
        )

    # Set tick label size
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize)

    ax.legend(fontsize=fontsize)
    ax.set_title(fig_title, fontsize=fontsize, y=1.03)
    ax.set_ylim([0, 0.85])
    return fig, ax
