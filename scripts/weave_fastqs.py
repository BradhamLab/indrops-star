import os
import itertools
import gzip

from Bio import SeqIO, SeqRecord

def combine_reads(r2, r4):
    """
    Combine R2 and R4 indrops v3 reads into barcode + UMI read.
    
    Parameters
    ----------
    r2 : Bio.SeqRecord
        Read containing the first half of the gel barcode. Should be from the R2
        fastq from an indrops V3 run.
    r4 : Bio.SeqRecord
        Read containing the second half of the gel barcode, the the UMI, and
        part of the poly-A tail. Should be from the R4 fastq from an indrops
        V3 run.
    
    Returns
    -------
    Bio.SeqRecord
        Interleafed barcode + UMI read where the first 16 characters are the
        cell barcode, and the next 6 are the UMI.
    """
    read_id = r2.id
    desc = r2.description.split(' ')[0] + ' bc+umi'
    # seq = r4.seq[8:8 + 6] + r2.seq + r4.seq[:8]
    # quality = {'phred_quality': r4.letter_annotations['phred_quality'][-6:] \
    #                           + r2.letter_annotations['phred_quality'] \
    #                           + r4.letter_annotations['phred_quality'][:-6]}
    seq = ''
    quality = {'phred_quality': []}
    for read in [r2, r4]:
        seq += read.seq
        quality['phred_quality'] += read.letter_annotations['phred_quality']
    return SeqRecord.SeqRecord(seq=seq, id=read_id, name=read_id,
                               description=desc, letter_annotations=quality)


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

def get_whitelist(barcode_file):
    """Get sequences from barcode file."""
    sequences = []
    with open(barcode_file, 'rU') as f:
        # iterate through each barcode (rstrip cleans string of whitespace)
        for line in f:
            barcode = line.rstrip()
            sequences.append(barcode)
    return sequences


def build_barcode_neighborhoods(sequences):
    """
    Create dictionary mapping observed barcodes to whitelist barcodes.

    Given a set of barcodes, produce sequences which can unambiguously be
    mapped to these barcodes, within 2 substitutions. If a sequence maps to 
    multiple barcodes, get rid of it. However, if a sequences maps to a bc1 with 
    1change and another with 2changes, keep the 1change mapping.

    Parameters
    ----------
        sequences : list
            List of sequences to generate neighborhoods for. 

    Returns
    -------
        dict (str, str)
            Dictionary where keys are mutated sequences and values corrected
            sequences. 

    References
    ----------
        Taken from original indrops repository: github.com/indrops/indrops
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


def get_library(seq, neighborhoods):
    """
    Return library name associated with given sequence.
    
    Parameters
    ----------
    seq : Bio.SeqRecord
        Library index read from R3 fastq of an indrops V3 run.
    neighborhoods : dict
        Dictionary of neighboring sequences for each library index. Where each
        key is a neighboring sequence, and each value is the associated library.
        name.
    
    Returns
    -------
    str
        Name of library associated with provided sequence. If sequence is not
        present, returns "ambig".
    """
    try:
        return neighborhoods[str(seq.seq)]
    except KeyError:
        return 'ambig'

    
def open_library_fastqs(libraries, prefix=''):
    """Open all fastq IO objects."""
    io_dict = {}
    for name in libraries.values():
        io_dict[name] = {'cdna': open(prefix + \
                                      "{}_cdna.fastq".format(name), 'w'),
                         'bc_umi': open(prefix + \
                                        "{}_bc_umi.fastq".format(name), 'w')}
        if name == 'ambig':
            io_dict[name]['index'] = open(prefix + \
                                          '{}_library_idx.fastq'.format(name),
                                          mode='w')
    # reads for non-mapped
    return io_dict


def close_fastqs(fastq_dict):
    """Close all fastq IO objects."""
    for each in fastq_dict:
        for every in fastq_dict[each].values():
            every.close()

            
def write_to_fastq(reads, fastqs):
    """
    Write reads to library-specific fastq files.

    Parameters
    ----------
        reads : dict
            Dictionary containing library-segregated reads.
        fastqs : list
            Dictionary holding library-specific fastq files. Output from
            `open_library_dict()`.
    """
    for library in fastqs.keys():
        SeqIO.write(reads[library]['cdna'], fastqs[library]['cdna'], 'fastq')
        SeqIO.write(reads[library]['bc_umi'], fastqs[library]['bc_umi'], 'fastq')
    if library == 'ambig':
            SeqIO.write(reads[library]['index'], fastqs[library]['index'],
                        'fastq')

            
def clear_reads(reads):
    """Clear current stored reads."""
    for library in reads:
        for key, value in reads[library].items():
            del value
            reads[library][key] = []

            
def parse_indrops_reads(r1_fastq, r2_fastq, r3_fastq, r4_fastq,
                        libraries, prefix='', nsubs=2):
    """
    Demultiplex Indrops V3 Reads.

    Parses the read files produced from Indrops V3 to 1.) segregate reads by
    library, and 2.) interleaf R2 and R4 files to produce a single read with 
    the complete cell barcode, UMI, and polyA tail.
    
    Parameters
    ----------
    r1_fastq : str
        File path to bio-read fastq of indrops v3 run.
    r2_fastq : str
        File path to fastq containing first half cell barcodes for an indrops
        v3 run.
    r3_fastq : str
        File path to fastq containing library indices for an indrops v3 run.
    r4_fastq : str
        File path to fastq containing second half of cell barcodes, UMI, and
        portion of the polyA tail from an indrops v3 run.
    libraries : dict
        Dictionary maping library index sequences to library names. Example
        `{"ATTAGAGG": "ASW-18hpf}
    prefix : str, optional
        Prefix/path to append to beginning of output files. Default is '', and
        no prefix is appended.
    nsubs : int, optional
        Allowable number of mismatches when matching libraries, by default 2.
    """
    if prefix != '' and not os.path.exists(os.path.basename(prefix)):
        os.makedirs(os.path.basename(prefix))
    if prefix != '' and prefix == os.path.basename(prefix):
        prefix = prefix + '/'
    # create IO objects to access reads
    fastq_handles = []
    reads = []
    for x in [r1_fastq, r2_fastq, r3_fastq, r4_fastq]:
        handle = gzip.open(x, 'rt')
        reads.append(SeqIO.parse(handle, 'fastq'))
        fastq_handles.append(handle)

    # construct dictionary mapping allowable library indices to library names
    library_neighborhoods = {}
    for index, name in libraries.items():
        library_neighborhoods.update(seq_neighborhood(index, name, nsubs))
    # create extra key '' for unmatched library reads
    libraries['ambig'] = "ambig"

    # open library segregated fastq files
    # dict[library] = {'<library>_cdna.fastq', -- bio read
    #                  '<library>_index.fastq', -- library index read
    #                  '<library>_bc_umi.fastq'} -- combined bc halves, umi, + polyA 
    fastqs = open_library_fastqs(libraries, prefix)
    library_dict = {name: {'cdna': [], 'index': [], 'bc_umi': []}\
                   for name in libraries.values()}
    # might be a way to use generators for speed -- unsure how to filter 
    # multiple read files using generator based off of r3, though.
    records = 0
    for r1, r2, r3, r4 in zip(*reads):
        # match library index to library name
        library = get_library(r3, library_neighborhoods)
        library_dict[library]['cdna'].append(r1)
        library_dict[library]['bc_umi'].append(combine_reads(r2, r4))
        if library == 'ambig':
            library_dict[library]['index'].append(r3)
        records += 1
        # write and clear reads every 100,000 records
        if records % 100000 == 0:
            # write current reads
            write_to_fastq(library_dict, fastqs)
            # clear written reads
            clear_reads(library_dict)
    # write any remaining reads
    write_to_fastq(library_dict, fastqs)
    # close all output file objects
    close_fastqs(fastqs)
    # close all input file objects
    for handle in fastq_handles:
        handle.close()


if __name__ == "__main__":
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        parse_indrops_reads(snakemake.input['r1'],
                            snakemake.input['r2'],
                            snakemake.input['r3'],
                            snakemake.input['r4'],
                            snakemake.params['libraries'],
                            prefix=snakemake.params['prefix'],
                            nsubs=snakemake.params['nsubs'])