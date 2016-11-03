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


def find_similar_residues(matrix=matinf.blosum62, threshold_min=0.5, threshold_max=0.0):
    """
    finds groups of similar residues based on a substitution matrix and score thresholds
    :param matrix: MatrixInfo variable containing common substitution matrices
    :param threshold_min: minimal score for inclusion in list of groups (exclusive)
    :param threshold_max: maximal score for inclusion in list of groups (inclusive)
    :return: groups and scores: [
                                    [
                                        ['R1', 'R2', ..., 'Rn'],        first group
                                        score
                                    ],
                                    [...],
                                    [
                                        ['R1', 'R2', ..., 'Rn'],        last group
                                        score
                                    ],
                                ]
    """
    no_max = False
    if threshold_max == 0.0:
        no_max = True

    letters = IUPACProtein.letters
    subs_matrix = matrix

    # reduce substitution matrix to residues in IUPACProtein.letters
    subs_matrix = {k: subs_matrix[k] for k in subs_matrix if k[0] in letters and k[1] in letters}

    # make binary correspondences
    residues_binary = {r: 2 ** letters.index(r) for r in letters}
    binary_residues = {str(2 ** letters.index(r)): r for r in letters}
    subs_matrix_binary = {str(residues_binary[m[0]] | residues_binary[m[1]]): subs_matrix[m] for m in subs_matrix}
    doubles = [2 ** i | 2 ** i for i in range(len(letters))]

    # pairs with negative values in substitution matrix
    neg_pairs = [residues_binary[m[0]] | residues_binary[m[1]] for m in subs_matrix if subs_matrix[m] < 0]
    neg_pairs = [p for p in neg_pairs if p not in doubles]

    # total combinations without doubles
    combs = [i for i in range(2 ** (len(letters)) - 1) if i not in doubles]

    # total positive combinations (all combinations not containing a negative or null pair)
    pos_combs = []
    for c in combs:
        t = True
        for n in neg_pairs:
            if c & n == n:
                t = False
                break
        if t:
            pos_combs.append(c)

    # combinations split up in single residue binaries
    pos_groups = [[b for b in doubles if b & r] for r in pos_combs]

    # positive groups and all possible combinations of residue pairs
    group_combs = [[g, [sum(c) for c in combinations(g, 2)]] for g in pos_groups if g]

    # positive groups with residues as 1 letter code, group score
    group_combs = [
        [
            [binary_residues[str(g)] for g in g[0]], sum([subs_matrix_binary[str(p)] for p in g[1]]) / len(g[1])
        ] for g in group_combs
    ]

    # positive groups with residues as 1 letter code, group score when group score
    # meets threshold requirements
    if no_max:
        group_combs_score = [[g[0], g[1]] for g in group_combs if g[1] > threshold_min]
    else:
        group_combs_score = [[g[0], g[1]] for g in group_combs if threshold_min < g[1] <= threshold_max]

    # positive groups with residues as 1 letter code, group score when group score
    # meets threshold requirements, only groups that are not a subset of another group
    final_groups = []
    for g, s in group_combs_score:
        tp = True
        for c in group_combs_score:
            if set(g).issubset(set(c[0])) and g != c[0]:
                tp = False
                break
        if tp:
            final_groups.append([g, s])
    return final_groups


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
