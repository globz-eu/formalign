from django.test import TestCase
import io
from Bio import SeqIO
from helper_funcs.helpers_bio import parse_fasta_alignment
from helper_funcs.helpers_test import file_to_string
from base.models import Alignment, Seqrecord

__author__ = 'Stefan Dieterle'


class SeqrecordModelTestCase(TestCase):
    """
    Tests Seqrecord model
    """

    def setUp(self):
        seq_input = io.StringIO(file_to_string('spa_align_clustal_omega.fasta'))
        data = SeqIO.parse(seq_input, 'fasta')
        self.seqs = [d.upper() for d in data]
        self.seq = Seqrecord.create(
                seq=self.seqs[0].seq,
                alphabet=str(self.seqs[0].alphabet)
        )

    def test_seqrecord_basic(self):
        """
        Tests the basic functionality of Seqrecord
        :return:
        """
        self.assertEqual(self.seq.alphabet, 'ExtendedIUPACProtein()')
        self.assertEqual(self.seq.seq, self.seqs[0].seq)


class AlignmentModelTestCase(TestCase):
    """
    Tests Alignment model
    """

    def setUp(self):
        name = 'A. tha. SPA family alignment'
        align_input = io.StringIO(file_to_string('spa_align_clustal_omega.fasta'))
        data = parse_fasta_alignment(align_input)
        self.alignment = Alignment.save_alignment_to_db(name, data)

    def test_alignment_basic(self):
        """
        Tests the basic functionality of Alignment
        :return:
        """
        self.assertEqual(self.alignment.name, 'A. tha. SPA family alignment')
        self.assertEqual(
                [seq.id for seq in self.alignment.seqs],
                ['AT1G53090.1', 'AT3G15354.1', 'AT2G46340.1', 'AT4G11110.1'])
