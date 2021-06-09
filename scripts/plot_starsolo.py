import re
import os

import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

import seaborn as sns

def log_split(line):
    """Split method for Log.final.out file produced by STAR."""
    return [x.strip() for x in line.strip().split(' |\t')]

def feature_split(line):
    """Split method for Feature.stats & Barcode.statsCa file produced by STARSolo."""
    return [x.strip() for x in  re.sub('\s+', '\t', line.strip()).split('\t')]

def stat_split(line):
    """Split method for Summary.csv file produced by STARSolo."""
    return [x.strip() for x in line.strip(',')]

def read_file(fn, split_method):
    """Parse STARSolo output files according to """
    log = {}
    with open(fn, 'r') as f:
        for line in f:
            splitter = split_method(line)
            if len(splitter) == 2:
                log[splitter[0]] = splitter[1].replace('%', '').strip()
    return log


def plot_star_log(logfiles, libraries, plotdir):
    # create two plots, one for n_reads for each library, one for alignment %s
    read_data = {}
    alignment_data = {}
    for (fn, lib) in zip(logfiles, libraries):
        log = read_file(fn, log_split)
        n_reads = int(log['Number of input reads'])
        unique = float(log['Uniquely mapped reads %'])
        multiple = float(log['% of reads mapped to multiple loci'])
        too_many = float(log['% of reads mapped to too many loci'])
        too_short = float(log['% of reads unmapped: too short'])
        mismatches = float(log['% of reads unmapped: too many mismatches'])
        unmapped = float(log['% of reads unmapped: other'])
        alignment_data[lib] = {'Unique': unique,
                               'Multiple': multiple + too_many, 
                               'Unmapped': too_short + unmapped + mismatches}

    align_df = pd.DataFrame.from_dict(alignment_data)
    read_df = pd.DataFrame.from_dict(read_data)

    align_ax = align_df.T.plot(kind='bar', stacked=True)
    align_ax.set_ylabel("% of Reads")
    align_ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    align_ax.set_title("Alignment Statistics")
    # plt.tight_layout()
    # plt.show()
    plt.savefig(os.path.join(plotdir, 'alignment.png'), bbox_inches='tight')


def plot_barcodes(statfiles, libraries, plotdir):
    barcode_data = {}
    umi_data = {}
    for (fn, lib) in zip(statfiles, libraries):
        stats = read_file(fn, feature_split)
        # print(stats)
        # ambiguous UMI + CB
        N_in_cb = int(stats['nNinCB'])
        N_in_UMI = int(stats['nNinUMI'])
        homopoly = int(stats['nUMIhomopolymer'])
        mismatch = int(stats['nMismatchToMultWL'])
        # no matches
        no_match = int(stats['nNoMatch'])

        matches = int(stats['nExactMatch'])
        barcode_data[lib] = {'Ambiguous': N_in_cb + mismatch,
                             'Matched': matches,
                             'No Match': no_match}
        umi_data[lib] = {'Homopolymers': homopoly, 'Ns present': N_in_UMI}

    bc_df = pd.DataFrame.from_dict(barcode_data)
    umi_df = pd.DataFrame.from_dict(umi_data)

    bc_ax = bc_df.T.plot(kind='bar', stacked=True)
    bc_ax.set_ylabel("# of Reads")
    bc_ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    bc_ax.set_title("Barcode Matching")
    plt.savefig(os.path.join(plotdir, 'barcodes.png'), bbox_inches='tight')


def additive_barplot(subset, ax=None):
    if ax is None:
        ax = plt.subplot()
    start = np.zeros(subset.shape[0])
    for col in subset.columns[::-1]:
        ax.bar(subset.index.values, subset[col].values - start,
               bottom=start, label=col)
        start = subset[col].values
    ax.legend(loc='upper left', bbox_to_anchor=(1,1))
    return ax


def plot_summaries(summaries, libraries, plotdir):
    dfs = []
    for fn, lib in zip(summaries, libraries):
        dfs.append(pd.read_csv(fn, index_col=0, header=None, names=[lib]).T)
    df = pd.concat(dfs)
    df['Valid Barcodes'] = df['Number of Reads'] * df['Reads With Valid Barcodes']
    df.rename(columns={"Number of Reads": 'Total'}, inplace=True)
    # df['Valid Barcodes'] = np.log10(df['barcodes'])  

    fig, axes = plt.subplots(figsize=(6, 9), nrows=3, ncols=1,
                             sharex=True, sharey=False)

    # n reads + valid barcodes
    ax = additive_barplot(df[['Total', 'Valid Barcodes']], ax=axes[0])
    ax.set_ylabel('# of Reads')
    ax.set_title('Reads')
    # saturation
    axes[1].bar(df.index, df['Sequencing Saturation'])
    axes[1].set_ylabel('Saturation')
    axes[1].set_title('Saturation')
    # mapping stats
    rename_cols = {x: x.replace("Reads Mapped to ", '')\
                   for x in df.columns if 'Reads Mapped' in x}
    df.rename(columns=rename_cols, inplace=True)
    ax = additive_barplot(df[rename_cols.values()], ax=axes[2])
    ax.set_title('Annotation Statistics')
    ax.set_ylabel('% of Reads')
    ax.set_xlabel('Library')
    plt.savefig(os.path.join(plotdir, 'summary.png'), bbox_inches='tight')

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        libraries = [x.split('/')[-2] for x in sorted(snakemake.input["star"])]
        if not os.path.exists(snakemake.params['plotdir']):
            os.makedirs(snakemake.params['plotdir'])
        plot_star_log(sorted(snakemake.input['star']),
                      libraries,
                      snakemake.params['plotdir'])
        plot_summaries(sorted(snakemake.input['summary']),
                       libraries,
                       snakemake.params['plotdir'])
