"""
=====================================================================
Formalign.eu format and display multiple sequence alignments
Copyright (C) 2016 Stefan Dieterle
e-mail: golgoths@yahoo.fr

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
=====================================================================
"""

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
