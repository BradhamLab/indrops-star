import os
import itertools

from Bio import SeqIO, SeqRecord

def combine_reads(r2, r4):
    """
    Combine R2 and R4 indrops v3 reads into barcode + UMI read.
    
    Parameters
    ----------
    r2 : Bio.SeqRecord
        Read containing the first half of the gel barcode. Should be from the R2
        fastq from an indrops V3 run.
    r4 : [type]
        Read containing the second half of the gel barcode, the the UMI, and
        part of the poly-A tail. Should be from the R4 fastq from an indrops
        V3 run.
    
    Returns
    -------
    Bio.SeqRecord
        Interleafed barcode + UMI read where the first 16 characters are the
        cell barcode, the next 6 are the UMI, and the remainder is part of the
        poly-A tail.
    """
    read_id = r2.id
    desc = r2.description.split(' ')[0] + ' R2R4'
    seq = ''
    quality = {'phred_quality': []}
    for read in [r2, r4]:
        seq += read.seq
        quality['phred_quality'] += read.letter_annotations['phred_quality']
    return SeqRecord.SeqRecord(seq=seq, id=read_id, name=read_id,
                               description=desc, letter_annotations=quality)


def seq_neighborhood(seq, name, nsubs=2):
    """
    Generate a sequence neighborhood for a given reference sequence

    Given a sequence, yield all sequences within `nsubs` substitutions of 
    that sequence by looping through each combination of base pairs within
    each combination of positions.

    Parameters
    ----------
    seq : str
        Reference sequence for neighbor generation. Example: library index.
    name : str
        Class/name for reference sequence. Example: library name.
    nsubs : int, optional
        Number of sub
    
    Returns
    -------
    dict
        Dictionary of neighborhood sequences such that each key is a neighboring
        sequence and each value is the provided name.

    References
    ----------
    Adapted from https://github.com/indrops/indrops
    """
    neighborhoods = {}
    for positions in itertools.combinations(range(len(seq)), nsubs):
    # yields all unique combinations of indices for nsubs mutations
        for subs in itertools.product(*("ATGCN",)*nsubs):
        # yields all combinations of possible nucleotides for strings of length
        # nsubs
            seq_copy = list(seq)
            for p, s in zip(positions, subs):
                seq_copy[p] = s
            neighborhoods[''.join(seq_copy)] = name
    return neighborhoods


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
    reads = [SeqIO.parse(x, format='fastq') for x in [r1_fastq, r2_fastq,
                                                      r3_fastq, r4_fastq]]
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
    # close all file objects
    close_fastqs(fastqs)


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