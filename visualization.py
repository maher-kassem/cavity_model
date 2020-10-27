import matplotlib.pyplot as plt
from typing import Tuple

def scatter_pred_vs_true(
    x,
    y,
    color: str,
    title: str,
    xlab: str = "Experimental Stability Change",
    ylab: str = "Predicted Stability Change",
    ylim: Tuple[float] = [-22., 22.],
    scatter_size: float = 0.5,
    fontsize: float = 15.0,
    labelpad: float = 10.0,
    title_gap: float = 1.03,
    figure_size: Tuple[float] = (5.0, 4.5),
    
):
    fig, ax = plt.subplots()
    ax.scatter(x, y, s=scatter_size, c=color)

    ax.set_ylim(ylim)
    max_abs_lim = max([abs(x) for x in ax.get_xlim()])
    ax.set_xlim([-max_abs_lim, max_abs_lim])

    ax.set_xlabel(xlab, fontsize=fontsize, labelpad=labelpad)
    ax.set_ylabel(ylab, fontsize=fontsize, labelpad=labelpad)
    
    ax.set_title(title, fontsize=fontsize, y=title_gap)
    
    fig.tight_layout()

    return fig, ax
