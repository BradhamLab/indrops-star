import re
import itertools

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt




def whitelist_barcodes(whitelist):
    """
    Calculate allowable barcodes from whitelist.

    Finds all 1-N substitutions that uniquely map to whitelist barcodes.

    Parameters
    ----------
        whitelist : str
            Path to file containing whitelist barcodes
    
    Returns
    -------
    dict
        Dictionary where keys are possible barcodes and values are whitelisted
        barcodes to map to.
    """
    barcodes = {}
    mapping = {}
    conflicting = {}
    # iterate through lines of whitelist barcodes
    with open(whitelist, 'r') as f:
        n_line = 0
        for line in f:
            n_line += 1
            print("{:0.2f}".format((n_line) / 147456), flush=True, end='\r') 
            bc = line.strip()
            # map barcode to barcode
            barcodes[bc] = bc
            # find sequences with a single substitution
            for i, base in itertools.product(range(len(bc)),
                                             ['A', 'T', 'C', 'G', 'N']):
                # check if substitution possible
                if bc[i] != base:
                    sub_bc = bc[:i] + base + bc[i + 1:]
                    if sub_bc not in conflicting and sub_bc not in mapping:
                        mapping[sub_bc] = bc
                    elif sub_bc in mapping:
                        conflicting[sub_bc] = ''
                        mapping.pop(sub_bc)
    # iterate through mapped barcodes, don't keep if substitution matches
    # whitelist barcode
    print('\n')
    print(f"mapping size: {len(mapping)}")
    print(f"conflicting substitutions: {len(conflicting)}")
    for k, v in mapping.items():
        if k not in barcodes:
            barcodes[k] = v
    return barcodes

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        whitelist = whitelist_barcodes(snakemake.input['whitelist'])
        bc_counts = {}
        count_re = re.compile('[0-9]*')
        bc_re = re.compile('[A-Z].*')
        with open(snakemake.input['counts'], 'r') as f:
            for line in f:
                line = line.strip()
                counts = int(count_re.search(line).group())
                bc = bc_re.search(line).group()
                if bc in whitelist:
                    mapped = whitelist[bc]
                    if mapped in bc_counts:
                        bc_counts[mapped] += counts
                    else:
                        bc_counts[mapped] = counts
        df = pd.DataFrame(bc_counts.items(),
                            columns=['barcode', 'reads'])
        df['whitelist'] = df.barcode.apply(lambda x: x in whitelist)
        df['log10(reads)'] = np.log10(df['reads'])
        df.to_csv(snakemake.output['csv'])
        subset = df[df['whitelist']]
        subset.describe().to_csv(snakemake.output['summary'])
        
        pdf = sns.distplot(subset['log10(reads)'])
        plt.savefig(snakemake.output['pdf'])
        plt.cla()
        plt.clf()
        plt.close()

        fig, ax = plt.subplots()
        # plot the cumulative histogram
        n, bins, patches = ax.hist(subset['log10(reads)'], 50, density=True,
                                   histtype='stepfilled',
                                   cumulative=True)
        # tidy up the figure
        ax.grid(True)
        # ax.legend(loc='right')
        ax.set_title('Read CDF')
        ax.set_xlabel('Log10(Reads)')
        ax.set_ylabel('Likelihood of occurrence')
        plt.savefig(snakemake.output['cdf'])  
                
