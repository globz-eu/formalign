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
import Bio.SubsMat.MatrixInfo as matinf
from Bio.Alphabet.IUPAC import IUPACProtein
from itertools import combinations

__author__ = 'Stefan Dieterle'


class AlignInfoSubsMat(AlignInfo.SummaryInfo):
    """
    subclasses SummaryInfo and overrides gap_consensus to use substitution matrices
    """


# def find_similar_residues(matrix=matinf.pam250, threshold_min=0.5, threshold_max=0.0, max_group_length=7):
#     """
#     finds groups of similar residues based on a substitution matrix and score thresholds
#     :param matrix: MatrixInfo variable containing common substitution matrices
#     :param threshold_min: minimal score for inclusion in list of groups (exclusive)
#     :param threshold_max: maximal score for inclusion in list of groups (inclusive)
#     :return: groups and scores: [
#                                     [
#                                         ['R1', 'R2', ..., 'Rn'],        first group
#                                         score
#                                     ],
#                                     [...],
#                                     [
#                                         ['R1', 'R2', ..., 'Rn'],        last group
#                                         score
#                                     ],
#                                 ]
#     """
#     no_max = False
#     if threshold_max == 0.0:
#         no_max = True
#     residues_binary = [(r, 2 ** IUPACProtein.letters.index(r)) for r in IUPACProtein.letters]
#     bins = [2 ** r for r in range(20)]
#     combs = [[r[0] for r in residues_binary if r[1] & i] for i in range(1, 1048576) if i not in bins]
#     combs = [set(c) for c in combs if 1 <= len(c) <= max_group_length]
#     combs.sort(key=len, reverse=False)
#     positive_scores = []
#     for c in combs:
#         if len(c) == 2:
#             try:
#                 score = matrix[tuple(c)]
#             except KeyError:
#                 score = matrix[tuple(c)[::-1]]
#             if score > 0:
#                 positive_scores.append([c, score])
#         else:
#             subset_scores = [sc for co, sc in positive_scores if c.issuperset(co)]
#             subsets = [co for co, sc in positive_scores if c.issuperset(co)]
#             if subset_scores:
#                 combinations_pairwise = [co for co in combinations(c, 2)]
#                 scores = [
#                     matrix[pair] for pair in matrix if (
#                         pair in combinations_pairwise or (pair[1], pair[0]) in combinations_pairwise
#                     )
#                 ]
#                 subset_scores.sort(reverse=True)
#                 if no_max:
#                     if threshold_min < (sum(scores) / len(scores)):
#                         score_comb = [c, sum(scores) / len(scores)]
#                         positive_scores = [p for p in positive_scores if p[0] not in subsets]
#                         positive_scores.append(score_comb)
#                 else:
#                     if threshold_min < (sum(scores) / len(scores)) <= threshold_max:
#                         score_comb = [c, sum(scores) / len(scores)]
#                         positive_scores = [p for p in positive_scores if p[0] not in subsets]
#                         positive_scores.append(score_comb)
#
#     positive_scores = [[list(co), sc] for co, sc in positive_scores if sc > threshold_min]
#     return positive_scores


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
