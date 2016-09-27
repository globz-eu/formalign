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

from django.test import TestCase

import io
from datetime import datetime, timedelta, timezone
from Bio import SeqIO
from Bio.Alphabet.IUPAC import ExtendedIUPACProtein
from Bio.Alphabet import Gapped

from helper_funcs.helpers_bio import parse_fasta_alignment
from helper_funcs.helpers_test import file_to_string

from base.models import Alignment, Seqrecord

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
        """
        Tests that an alignment can be retrieved directly from the database after having been saved to the database
        """
        save = Alignment.objects.create_alignment(self.name, self.data)
        alignment = Alignment.objects.all()
        self.assertEqual(alignment[0].name, 'A. tha. SPA family alignment', alignment[0].name)
        self.assertEqual(
            [('NP_175717', 0), ('NP_683567', 1), ('NP_182157', 2), ('NP_192849', 3)],
            [(seq.seq_id, seq.display_order) for seq in alignment[0].seqs.all()],
            [(seq.seq_id, seq.display_order) for seq in alignment[0].seqs.all()]
        )
        self.assertEqual(save.id, alignment[0].pk, save.id)

    def test_save_alignment_adds_created_timestamp(self):
        """
        Tests that when an alignment is saved to the database it contains a created time stamp
        """
        Alignment.objects.create_alignment(self.name, self.data)
        alignment = Alignment.objects.all()
        self.assertTrue(datetime.now(timezone.utc) - alignment[0].created < timedelta(minutes=1))

    def test_get_multipleseqalignment_object_from_db(self):
        """
        Tests that an alignment can be retrieved as a MultipleSeqAlignment object after having been saved to the
        database
        """
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
            ['NP_175717', 'NP_683567', 'NP_182157', 'NP_192849'],
            [seq.id for seq in mulseqal],
            [seq.id for seq in mulseqal],
        )

    def test_save_and_get_sanity_check(self):
        """
        Tests that an alignment can be retrieved intact after having been saved to the database
        :return:
        """
        align_pk = Alignment.objects.create_alignment(self.name, self.data).id
        mulseqal = Alignment.objects.get_alignment(align_pk)
        self.assertEqual(
            ['NP_175717', 'NP_683567', 'NP_182157', 'NP_192849'],
            [seq.id for seq in mulseqal],
            [seq.id for seq in mulseqal],
        )
        self.assertEqual(str(mulseqal), str(self.data), str(mulseqal))
