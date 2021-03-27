from collections import defaultdict
from itertools import combinations, product

# ---------------------- Sequence Neighborhood Matching ------------------------
def seq_neighborhood(seq, n_subs=1):
    """
    Generate all n-substitution strings from sequence.

    Given a sequence, yield all sequences within n_subs substitutions of 
    that sequence by looping through each combination of base pairs within
    each combination of positions.

    Parameters
    ----------
        seq : str
            Base sequence to mutate over.
        n_subs : int
            Number of allowable substitutions.
    
    Returns
    -------
        Generator
            All possible n-substitution strings
    
    References
    ----------
        Taken from original indrops repository: github.com/indrops/indrops
    """
    for positions in combinations(range(len(seq)), n_subs):
    # yields all unique combinations of indices for n_subs mutations
        for subs in product(*("ATGCN",)*n_subs):
        # yields all combinations of possible nucleotides for strings of length
        # n_subs
            seq_copy = list(seq)
            for p, s in zip(positions, subs):
                seq_copy[p] = s
            yield ''.join(seq_copy)


def build_sequence_neighborhoods(sequences):
    """
    Create dictionary mapping allowable mutated sequence to expected sequence.

    Given a set of sequences, produce mutated sequences which can unambiguously
    be mapped to expected sequences, within 2 substitutions. If a mutatation
    maps to  multiple sequences, get rid of it. However, if a mutation maps to
    an expected sequence with 1change and another with 2changes,
    keep the 1change mapping.

    Parameters
    ----------
        sequences : list, generator
            List of sequences to generate neighborhoods for. 

    Returns
    -------
        dict (str, str)
            Dictionary where keys are mutated sequences and values corrected
            sequences. 

    References
    ----------
        Modified from original indrops repository: github.com/indrops/indrops
    """

    # contains all mutants that map uniquely to a barcode
    clean_mapping = dict()

    # contain single or double mutants 
    mapping1 = defaultdict(set)
    mapping2 = defaultdict(set)
    
    #Build the full neighborhood and iterate through barcodes
    for seq in sequences:
        # each barcode obviously maps to itself uniquely
        clean_mapping[seq] = seq

        # for each possible mutated form of a given barcode, either add
        # the origin barcode into the set corresponding to that mutant or 
        # create a new entry for a mutant not already in mapping1
        # eg: barcodes CATG and CCTG would be in the set for mutant CTTG
        # but only barcode CATG could generate mutant CANG
        for n in seq_neighborhood(seq, 1):
            mapping1[n].add(seq)
        
        # same as above but with double mutants
        for n in seq_neighborhood(seq, 2):
            mapping2[n].add(seq)   
    
    # take all single-mutants and find those that could only have come from one
    # specific barcode
    for k, v in mapping1.items():
        if k not in clean_mapping:
            if len(v) == 1:
                clean_mapping[k] = list(v)[0]
    
    for k, v in mapping2.items():
        if k not in clean_mapping:
            if len(v) == 1:
                clean_mapping[k] = list(v)[0]
    del mapping1
    del mapping2
    return clean_mapping