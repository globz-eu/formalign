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

import io
import os
import pickle

from unittest import TestCase as UniTestCase
from django.test import TestCase

from helper_funcs.helpers_bio import parse_fasta_alignment, consensus_add
# from helper_funcs.helpers_bio import find_similar_residues
from helper_funcs.helpers_test import file_to_string
from helper_funcs.helpers_format import split_lines, split_lines_in_blocks, annotate

from Bio.Alphabet.IUPAC import ExtendedIUPACProtein
from Bio.Alphabet import Gapped
from Bio.Align import AlignInfo
# import Bio.SubsMat.MatrixInfo as matinf

__author__ = 'Stefan Dieterle'


class FormatHelpersTestCase(TestCase):
    """
    Tests for seq-display and align_display views helper functions
    """

    def setUp(self):
        """
        Creates an alignment from ser_thr_kin_short in the db
        """
        align_input = io.StringIO(file_to_string('ser_thr_kin_short.fasta'))
        data = parse_fasta_alignment(align_input)
        for d in data:
            d.seq.alphabet = Gapped(ExtendedIUPACProtein())
        self.alignment = data
        align_input_a = io.StringIO(file_to_string('protein_annotate_test.fasta'))
        data_a = parse_fasta_alignment(align_input_a)
        for d in data_a:
            d.seq.alphabet = Gapped(ExtendedIUPACProtein())
        self.alignment_a = data_a

    def test_consensus_get_returns_consensus_to_alignment(self):
        """
        Tests that consensus_add function returns correct default consensus
        """
        alignment = self.alignment
        cons_seq = AlignInfo.SummaryInfo(alignment).gap_consensus()
        cons_got = consensus_add(alignment)[-1]
        self.assertEqual(cons_seq, cons_got.seq, cons_got.seq)
        self.assertEqual('consensus 70%', cons_got.id, cons_got.id)

    def test_split_lines_alignment_splits_blocks_of_80_chars(self):
        """
        Tests that split_lines function splits blocks of alignment sequences in lines of 80 characters by default
        """
        alignment = self.alignment
        split_alignment = split_lines(alignment, split_type='alignment')
        for seq in split_alignment[:-1]:
            self.assertEqual([80 for i in range(len(seq))], [len(s) for s in seq])
        for seq in split_alignment[-1]:
            self.assertTrue(len(seq) <= 80)

    def test_split_lines_sequence_splits_lines_of_80_chars(self):
        """
        Tests that split_lines function splits sequences in lines of 80 characters by default
        """
        alignment = self.alignment
        split_sequence = split_lines(alignment, split_type='sequence')
        for seq in split_sequence:
            self.assertEqual([80 for i in range(len(seq['seq']) - 1)], [len(s) for s in seq['seq'][:-1]])
            self.assertTrue(len(seq['seq'][-1]) <= 80)

    def test_split_lines_in_blocks_returns_correct_object(self):
        """
        Tests that lines are split in blocks of 3 characters, tests that split_lines_in_blocks returns the correct
        object
        """
        split_seq_file = os.path.join(
            os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'test_data', 'split_seq.pckl'
        )
        split_sequence = pickle.load(open(split_seq_file, 'rb'))
        split_blocks = split_lines_in_blocks(split_sequence, block_size=3)
        split_blocks_li_exp = [
            [
                ['seq1', [[(0, 'C'), (1, 'N'), (0, 'G')], [(1, 'R'), (0, 'S'), (1, 'S')]]],
                ['seq2', [[(0, 'E'), (1, 'N'), (0, 'G')], [(1, 'R'), (0, 'G'), (1, 'S')]]],
                ['seq3', [[(0, 'D'), (0, 'R'), (0, 'D')], [(1, 'R'), (0, 'C'), (1, 'S')]]],
                ['seq4', [[(0, 'R'), (1, 'N'), (0, '-')], [(0, '-'), (0, '-'), (1, 'S')]]],
                ['consensus 70%', [[(0, 'X'), (1, 'N'), (0, 'X')], [(1, 'R'), (0, 'X'), (1, 'S')]]]
            ],
            [
                ['seq1', [[(0, '-'), (0, '-'), (0, '-')], [(0, '-')]]],
                ['seq2', [[(0, 'M'), (0, 'I'), (0, 'V')], [(0, 'S')]]],
                ['seq3', [[(0, '-'), (0, '-'), (0, '-')], [(0, '-')]]],
                ['seq4', [[(0, '-'), (0, '-'), (0, '-')], [(0, '-')]]],
                ['consensus 70%', [[(0, '-'), (0, '-'), (0, '-')], [(0, '-')]]]
            ]
        ]
        split_blocks_li = [
            [
                [line[0], [[(annot, res) for (annot, res) in res_block] for res_block in line[1]]
                 ] for line in seq_block
                ] for seq_block in split_blocks
            ]
        for seq_block in split_blocks_li:
            for line in seq_block:
                res_blocks = line[1]
                for res_block in res_blocks:
                    length = len([(annot, res) for (annot, res) in res_block])
                    self.assertTrue(length <= 3, length)
        self.assertEqual(split_blocks_li_exp, split_blocks_li, split_blocks_li)

    def test_annotate_returns_correct_default_consensus_annotation(self):
        """
        Tests that consensus is correctly annotated for display
        """
        exp_annot = [0, 1, 0, 1, 0, 1, 0, 0, 0, 0]
        alignment = consensus_add(self.alignment_a)
        annot_alignment = annotate(alignment)
        self.assertEqual(
            exp_annot,
            annot_alignment[-1].letter_annotations['eq'],
            annot_alignment[-1].letter_annotations['eq']
        )

    def test_annotate_returns_correct_sequence_annotations(self):
        """
        Tests that sequences are correctly annotated for display
        """
        exp_annot = [
            [0, 1, 0, 1, 0, 1, 0, 0, 0, 0],
            [0, 1, 0, 1, 0, 1, 0, 0, 0, 0],
            [0, 0, 0, 1, 0, 1, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 1, 0, 0, 0, 0],
            [0, 1, 0, 1, 0, 1, 0, 0, 0, 0]
        ]
        alignment = consensus_add(self.alignment_a)
        annot_alignment = annotate(alignment)
        for i, e in enumerate(exp_annot):
            self.assertEqual(
                e,
                annot_alignment[i].letter_annotations['eq'],
                annot_alignment[i].letter_annotations['eq']
            )


class BioHelpers(UniTestCase):
    """
    tests for parse_fasta_alignment and consensus_add functions
    """

    def test_parse_fasta_alignment_returns_expected_object(self):
        """
        tests that parse_fasta_alignment returns the expected object
        """
        align = io.StringIO(file_to_string('ser_thr_kin_short.fasta'))
        parsed = parse_fasta_alignment(align)
        self.assertEqual(
            ['DMD401_1-640', 'CER09D1_11-435', 'EGFR', 'DMDPR2_1-384'],
            [p.description for p in parsed],
            [p.description for p in parsed]
        )

    # def test_find_similar_residues_returns_the_correct_object(self):
    #     """
    #     tests that find_similar_residues returns the correct object
    #     """
    #     blosum_62_high = find_similar_residues(matinf.blosum62, threshold_min=0.6, max_group_length=6)
    #     blosum_62_high = [[set(co), sc] for co, sc in blosum_62_high]
    #     self.assertEqual(
    #         [[{'H', 'N'}, 1], [{'A', 'T', 'S'}, 0.6666666666666666], [{'E', 'K', 'Q', 'D'}, 0.8333333333333334],
    #          [{'E', 'N', 'Q', 'D'}, 0.8333333333333334], [{'E', 'N', 'D', 'S'}, 0.6666666666666666],
    #          [{'W', 'Y', 'F', 'H'}, 0.8333333333333334], [{'E', 'K', 'N', 'Q', 'R'}, 0.7],
    #          [{'M', 'L', 'V', 'A', 'I'}, 0.7], [{'M', 'L', 'F', 'I', 'V'}, 0.9]],
    #         blosum_62_high, blosum_62_high)
    #     pam_250_high = find_similar_residues(matinf.pam250, threshold_min=1.3, max_group_length=7)
    #     pam_250_high = [[set(co), sc] for co, sc in pam_250_high]
    #     self.assertEqual(
    #         [[{'R', 'W'}, 2], [{'G', 'D', 'E'}, 1.3333333333333333], [{'Y', 'H', 'F'}, 1.6666666666666667],
    #          [{'Q', 'D', 'H', 'R'}, 1.3333333333333333], [{'N', 'Q', 'R', 'K', 'H'}, 1.4],
    #          [{'M', 'I', 'L', 'F', 'V'}, 1.8], [{'H', 'N', 'Q', 'D', 'K', 'E'}, 1.3333333333333333]],
    #         pam_250_high, pam_250_high)
