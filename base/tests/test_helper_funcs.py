import io
import os
import pickle
from unittest import TestCase as UniTestCase

import Bio.SubsMat.MatrixInfo as matinf
from Bio.Align import AlignInfo
from Bio.Alphabet import Gapped
from Bio.Alphabet.IUPAC import ExtendedIUPACProtein
from django.test import TestCase
from helper_funcs.bio.helpers import find_similar_residues
from helper_funcs.bio.helpers import parse_fasta_alignment, consensus_add
from helper_funcs.helpers_format import split_lines, split_lines_in_blocks, annotate
from helper_funcs.helpers_test import file_to_string


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

    def test_find_similar_residues_returns_the_correct_object(self):
        """
        tests that find_similar_residues returns the correct object (groups and scores)
        """
        blosum_62_high = find_similar_residues(threshold_min=0.5)
        self.assertEqual(
            [[['F', 'I', 'L', 'M'], 0.8333333333333334], [['H', 'N'], 1.0], [['E', 'H', 'Q'], 0.6666666666666666],
             [['E', 'K', 'N', 'Q', 'R'], 0.7], [['E', 'K', 'Q', 'S'], 0.6666666666666666],
             [['D', 'E', 'N', 'Q', 'S'], 0.6], [['A', 'S', 'T'], 0.6666666666666666],
             [['N', 'S', 'T'], 0.6666666666666666], [['I', 'L', 'M', 'V'], 1.6666666666666667], [['H', 'Y'], 2.0],
             [['F', 'W', 'Y'], 2.0]],
            blosum_62_high, blosum_62_high)
        blosum_62_low = find_similar_residues(matinf.blosum62, threshold_min=0, threshold_max=0.5)
        self.assertEqual(
            [[['F', 'I', 'M'], 0.3333333333333333], [['E', 'K', 'N', 'R'], 0.5], [['E', 'H', 'N', 'Q', 'R'], 0.4],
             [['A', 'G', 'S'], 0.3333333333333333], [['G', 'N', 'S'], 0.3333333333333333],
             [['D', 'N', 'Q', 'S'], 0.3333333333333333], [['E', 'K', 'N', 'Q', 'S'], 0.5]],
            blosum_62_low, blosum_62_low)
        pam_250_high = find_similar_residues(matinf.pam250, threshold_min=0.5)
        self.assertEqual(
            [[['F', 'I', 'L', 'M'], 1.8333333333333333], [['A', 'D', 'E', 'N', 'Q'], 1.1],
             [['D', 'E', 'H', 'K', 'N', 'Q'], 1.3333333333333333], [['K', 'M', 'R'], 1.0],
             [['H', 'K', 'N', 'Q', 'R'], 1.4], [['H', 'N', 'P', 'Q', 'R'], 0.9],
             [['K', 'N', 'R', 'S'], 0.8333333333333334], [['A', 'D', 'E', 'G', 'N', 'S', 'T'], 0.6190476190476191],
             [['D', 'E', 'K', 'N', 'S', 'T'], 0.6], [['A', 'G', 'N', 'P', 'S', 'T'], 0.5333333333333333],
             [['I', 'L', 'M', 'V'], 2.6666666666666665], [['I', 'T', 'V'], 1.3333333333333333], [['R', 'W'], 2.0],
             [['F', 'W', 'Y'], 2.3333333333333335]],
            pam_250_high, pam_250_high)
