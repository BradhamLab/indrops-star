from Bio import SeqIO, SeqRecord
import os
import itertools

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


def seq_neighborhood(seq, name, n_subs=2):
    """
    Generate a sequence neighborhood for a given reference sequence

    Given a sequence, yield all sequences within `n_subs` substitutions of 
    that sequence by looping through each combination of base pairs within
    each combination of positions.

    Parameters
    ----------
    seq : str
        Reference sequence for neighbor generation. Example: library index.
    name : str
        Class/name for reference sequence. Example: library name.
    n_subs : int, optional
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
    for positions in itertools.combinations(range(len(seq)), n_subs):
    # yields all unique combinations of indices for n_subs mutations
        for subs in itertools.product(*("ATGCN",)*n_subs):
        # yields all combinations of possible nucleotides for strings of length
        # n_subs
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

def open_library_fastqs(libraries, out_dir=''):
    """Open all fastq IO objects."""
    io_dict = {}
    for name in libraries.values():
        io_dict[name] = {'cdna': open(os.path.join(out_dir,
                                      "{}_cdna.fastq".format(name)), 'w'),
                         'bc_umi': open(os.path.join(out_dir,
                                        "{}_bc_umi.fastq".format(name)), 'w')}
        if name == 'ambig':
            io_dict[name]['index'] = open(os.path.join(out_dir,
                                          '{}_library_idx.fastq'.format(name)),
                                          mode='w')
    # reads for non-mapped
    return io_dict

def close_fastqs(fastq_dict):
    """Close all fastq IO objects."""
    for each in fastq_dict:
        for every in fastq_dict[each].values():
            every.close()

def parse_indrops_reads(r1_fastq, r2_fastq, r3_fastq, r4_fastq,
                        libraries, out_dir='', n_subs=2):
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
    out_dir : str, optional
        Where to write weaved fastq files. Current directory by default.
    n_subs : int, optional
        Allowable number of mismatches when matching libraries, by default 2.
    """
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    # create IO objects to access reads 
    reads = [SeqIO.parse(x, format='fastq') for x in [r1_fastq, r2_fastq,
                                                      r3_fastq, r4_fastq]]
    # construct dictionary mapping allowable library indices to library names
    library_neighborhoods = {}
    for index, name in libraries.items():
        library_neighborhoods.update(seq_neighborhood(index, name, n_subs))
    # create extra key '' for unmatched library reads
    libraries['ambig'] = "ambig"

    # open library segregated fastq files
    # --> dict[library] = {'<library>_cdna.fastq', '<library>_bc_umi.fastq'}
    fastqs = open_library_fastqs(libraries, out_dir)
    # maybe barcode correct just in case?
    for r1, r2, r3, r4 in zip(*reads):
        # match library index to library name
        library = get_library(r3, library_neighborhoods)
        # combine barcode halves + umi + polyA
        bc_umi = combine_reads(r2, r4)
        # write bio read to library segregated fastq
        SeqIO.write(r1, fastqs[library]['cdna'], 'fastq')
        # write bc+umi read to library segregated fastq
        SeqIO.write(bc_umi, fastqs[library]['bc_umi'], 'fastq')
        if library == 'ambig':
            SeqIO.write(r3, fastqs[library]['index'], 'fastq')

    # close all file objects
    close_fastqs(fastqs)