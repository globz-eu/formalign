from Bio import SeqIO

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
