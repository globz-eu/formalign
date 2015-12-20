from django.test import TestCase

import io

from Bio import SeqIO
from Bio.Alphabet.IUPAC import ExtendedIUPACProtein
from Bio.Alphabet import Gapped

from helper_funcs.helpers_bio import parse_fasta_alignment
from helper_funcs.helpers_test import file_to_string

from base.models import Alignment, Seqrecord  #, save_alignment_to_db, get_multipleseqalignment_object_from_db

__author__ = 'Stefan Dieterle'


class SeqrecordModelTestCase(TestCase):
    """
    Tests Seqrecord model
    """

    def setUp(self):
        seq_input = io.StringIO(file_to_string('spa_protein_alignment.fasta'))
        data = SeqIO.parse(seq_input, 'fasta')
        self.seqs = [d.upper() for d in data]
        for s in self.seqs:
            s.seq.alphabet = ExtendedIUPACProtein()
        self.seq = Seqrecord.objects.create(
                seq=str(self.seqs[0].seq),
                alphabet=str(self.seqs[0].seq.alphabet),
                seq_id=str(self.seqs[0].id),
                name=str(self.seqs[0].name),
                description=self.seqs[0].description,
        )

    def test_seqrecord_basic(self):
        """
        Tests the basic functionality of Seqrecord
        :return:
        """
        self.assertEqual(self.seq.alphabet, 'ExtendedIUPACProtein()', self.seq.alphabet)
        self.assertEqual(self.seq.seq, str(self.seqs[0].seq), self.seq.seq)


class AlignmentModelTestCase(TestCase):
    """
    Tests Alignment model
    """

    def setUp(self):
        self.name = 'A. tha. SPA family alignment'
        align_input = io.StringIO(file_to_string('spa_protein_alignment.fasta'))
        self.data = parse_fasta_alignment(align_input)
        alphabet = Gapped(ExtendedIUPACProtein())
        for a in self.data:
            a.seq.alphabet = alphabet
        self.data._alphabet = alphabet

    def test_alignment_basic(self):
        """
        Tests the basic functionality of Alignment
        :return:
        """
        alignment = Alignment.objects.create(name=self.name)
        for s in self.data:
            seqrec = Seqrecord.objects.create(
                    seq=str(s.seq),
                    alphabet=str(s.seq.alphabet),
                    seq_id=s.id,
                    name=s.name,
                    description=s.description
            )
            alignment.seqs.add(seqrec)
        self.assertEqual(alignment.name, 'A. tha. SPA family alignment')
        self.assertEqual(
                [seq.seq_id for seq in alignment.seqs.all()],
                ['NP_175717', 'NP_683567', 'NP_182157', 'NP_192849']
        )

    def test_save_alignment_to_db(self):
        save = Alignment.objects.create_alignment(self.name, self.data)
        alignment = Alignment.objects.all()
        self.assertEqual(alignment[0].name, 'A. tha. SPA family alignment', alignment[0].name)
        self.assertEqual(
                [seq.seq_id for seq in alignment[0].seqs.all()],
                ['NP_175717', 'NP_683567', 'NP_182157', 'NP_192849']
        )
        self.assertEqual(save.id, alignment[0].pk, save.id)

    def test_get_multipleseqalignment_object_from_db(self):
        alignment = Alignment.objects.create(name=self.name)
        for s in self.data:
            seqrec = Seqrecord.objects.create(
                    seq=str(s.seq),
                    alphabet="Gapped(ExtendedIUPACProtein(), '-')",
                    seq_id=s.id,
                    name=s.name,
                    description=s.description
            )
            alignment.seqs.add(seqrec)
        mulseqal = Alignment.objects.get_alignment(alignment.pk)
        self.assertEqual(
                [seq.id for seq in mulseqal],
                ['NP_175717', 'NP_683567', 'NP_182157', 'NP_192849']
        )

    def test_save_and_get_sanity_check(self):
        align_pk = Alignment.objects.create_alignment(self.name, self.data).id
        mulseqal = Alignment.objects.get_alignment(align_pk)
        self.assertEqual(
                [seq.id for seq in mulseqal],
                ['NP_175717', 'NP_683567', 'NP_182157', 'NP_192849']
        )
        self.assertEqual(str(mulseqal), str(self.data), str(mulseqal))
