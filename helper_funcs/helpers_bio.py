from Bio import SeqIO, AlignIO
from Bio.Align import AlignInfo
from io import StringIO

__author__ = 'Stefan Dieterle'


def parse_fasta(fasta):
        """
        Parses fasta file or file-like object into list of dicts of metadata and sequences using SeqIO.parse
        :param fasta: fasta file or file-like object
        :return: fasta_list = [{'meta': 'sequence meta', 'seq': 'SEQUENCE'} ... ]
        """
        fasta_parse = SeqIO.parse(fasta, 'fasta')
        fasta_dict_list = [
            {
               'meta': f.description,
               'seq':  str(f.seq).upper()
            } for f in fasta_parse
            ]
        return fasta_dict_list


def parse_fasta_alignment(fasta):
    """
    Parses fasta file or file-like object into list of dicts of metadata and sequences using AlignIO.read()
    :param fasta: fasta file or file-like object
    :return: fasta_list = [{'meta': 'sequence meta', 'seq': 'SEQUENCE'} ... ]
    """
    fasta_parse = AlignIO.read(fasta, 'fasta')
    for f in fasta_parse:
        f.seq = str(f.seq).upper()
    return fasta_parse


if __name__ == '__main__':
    from helper_funcs.helpers_test import file_to_string
    fasta_string = file_to_string('spa_align_clustal_omega.fasta')
    fasta = StringIO(fasta_string)
    parsed = parse_fasta_alignment(fasta)
    fasta = StringIO(fasta_string)
    alio_parsed = AlignIO.read(fasta, 'fasta')
    for p in parsed:
        print('seq: ', p.seq)
        print('meta: ', p.description)
    for a in alio_parsed:
        print('record: ', a)
