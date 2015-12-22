from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment, AlignInfo
from Bio.SeqRecord import SeqRecord

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


def consensus_add(alignment):
    """
    Calculates the appropriate consensus sequence for the given alignment
    :param alignment: MultipleSeqAlignment object
    :return: alignment with appended consensus as a SeqRecord object
    """
    cons_seqrec = SeqRecord(
                AlignInfo.SummaryInfo(alignment).gap_consensus(),
                id='consensus 70%',
                name='consensus 70%',
                description='consensus 70%',
        )
    alignment.append(cons_seqrec)
    return alignment
