from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

__author__ = 'Stefan Dieterle'


def parse_fasta_alignment(fasta):
    """
    Parses fasta file or file-like object into list of dicts of metadata and sequences using AlignIO.read()
    :param fasta: fasta file or file-like object
    :return: fasta_parsed = MultipleSeqAlignment object
    """
    fasta_parsed = AlignIO.read(fasta, 'fasta')
    fasta_parsed = MultipleSeqAlignment([f.upper() for f in fasta_parsed], annotations=fasta_parsed.annotations)
    return fasta_parsed
