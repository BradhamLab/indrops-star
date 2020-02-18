
import subprocess as sbp
import os
import itertools
import numpy as np

# function to get genomeChrBinNBits parameter for STAR alignment.
def estimate_STAR_ChrBinNbits(genome_file, read_length):
    """
    Estimate the `ChrBinNBits` parameter for genome indexing in STAR
    Estimate the `ChrBinNBits` parameter for genome indexing in STAR. Value
    must be estimated due to memory constraints caused by the large number
    of scaffolds present in some genomes (i.e. the LV genome). If estimation
    is unnecessary, flag `star_est_ChrBinNbits: False` in configuration file.
    Args:
        genome_file (string): path to fasta file containing genome reference
            sequences.
        read_length (int): length of reads from RNAseq experiment.
    Return:
        (int) new value for scaling RAM consumption
    References:
    https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf (p. 7)
    https://github.com/alexdobin/STAR/issues/103
    """
    len_call = 'grep -v ">" {} | wc | awk '.format(genome_file)\
               + "'{print $3-$1}'"
    n_ref_call = 'grep "^>" {} | wc -l'.format(genome_file)

    return_values = [None, None]
    for i, call in enumerate([len_call, n_ref_call]):
        p = sbp.Popen(call, stdin=sbp.PIPE, stdout=sbp.PIPE, stderr=sbp.PIPE,
                      shell=True)
        output, err = p.communicate()
        if p.returncode == 0:
            return_values[i] = int(output.strip())
        else:
            raise OSError(err)
    estimate = max([int(np.log2(return_values[0] / return_values[1])),
                    int(np.log2(read_length))])
    return min(18, estimate)