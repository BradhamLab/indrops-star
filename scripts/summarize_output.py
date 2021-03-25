import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

def sorted_umis_plot(umis, ax=None):
    if ax is None:
        ax = plt.subplot()
    ax.plot(np.arange(len(umis)) / len(umis) * 100, umis, color='black', linewidth=1.5)
    ax.set_ylabel('Counts')
    ax.set_xlabel('UMIs (percent)')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    return ax

def distplot(umis, ax=None):
    if ax is None:
        ax = plt.subplot()
    counts = ax.hist(umis, bins=30, histtype='stepfilled', color='lightgray',
                      edgecolor='black', linewidth=1.5, density=False)
    ax.set_xlabel('Counts')
    ax.set_ylabel('UMIs')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    return ax

def cumulative_plot(umis, ax=None):
    if ax is None:
        ax = plt.subplot()
    cumulative = ax.hist(umis, cumulative=True, bins=100, density=True,
                         histtype='stepfilled', color='lightgray', edgecolor='black',
                         linewidth=1.5)
    ax.set_xlabel('Counts')
    ax.set_ylabel('UMIs (percent)')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    return ax

def plot_umis(umis, library=''):
    fig = plt.figure()
    grid = plt.GridSpec(2, 2, hspace=0.5, wspace=0.5)
    rank_plot = sorted_umis_plot(umis, fig.add_subplot(grid[0, 0]))
    hist = distplot(umis, fig.add_subplot(grid[0, 1]))
    cumulative = cumulative_plot(umis, fig.add_subplot(grid[1, :]))
    title = 'UMI Distribution'
    if library != '':
        title = f'{library} {title}'
    plt.suptitle(title, x=0, horizontalalignment='left')
    return fig