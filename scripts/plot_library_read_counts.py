from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os

if __name__ == "__main__":
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        head_dir = os.path.split(os.path.dirname(__file__))[0]
        plt.style.use(os.path.join(head_dir, 'files', 'figures.mplstyle'))
        read_counts = pd.read_csv(snakemake.input['csv'],
                                  names=['Read.Counts', 'Library'])
        read_counts["Log10(Read Counts)"] = np.log10(read_counts['Read.Counts'])
        plt.subplots(figsize=(13, 8))
        sns.barplot(x='Library', y='Log10(Read Counts)', data=read_counts)
        plt.tight_layout()
        plt.savefig(snakemake.output['img'])
        plt.close()