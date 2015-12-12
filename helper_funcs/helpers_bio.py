from Bio import AlignIO

__author__ = 'Stefan Dieterle'


def parse_fasta_alignment(fasta):
    """
    Parses fasta file or file-like object into list of dicts of metadata and sequences using AlignIO.read()
    :param fasta: fasta file or file-like object
    :return: fasta_parse = AlignIO object
    """
    fasta_parse = AlignIO.read(fasta, 'fasta')
    for f in fasta_parse:
        f.seq = str(f.seq).upper()
    return fasta_parse
