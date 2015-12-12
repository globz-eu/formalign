from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

__author__ = 'Stefan Dieterle'


def parse_fasta_alignment(fasta):
    """
    Parses fasta file or file-like object into list of dicts of metadata and sequences using AlignIO.read()
    :param fasta: fasta file or file-like object
    :return: fasta_parse = AlignIO object
    """
    fasta_parsed = AlignIO.read(fasta, 'fasta')
    fasta_parsed = MultipleSeqAlignment([f.upper() for f in fasta_parsed], annotations=fasta_parsed.annotations)
    return fasta_parsed


if __name__ == '__main__':
    # from helper_funcs.helpers_test import file_to_string
    # fasta = file_to_string('spa_align_clustal_omega.fasta')
    al = parse_fasta_alignment('../test_data/spa_align_clustal_omega.fasta')
    # for a in al:
    #     print(a.upper())
    print(al)
