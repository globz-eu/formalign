import io
import re
from datetime import datetime, timedelta, timezone

from Bio import SeqIO
from Bio.Alphabet import Gapped
from Bio.Alphabet.IUPAC import ExtendedIUPACProtein
from base.models import Alignment, Seqrecord
from django.test import TestCase
from helper_funcs.bio.helpers import parse_fasta_alignment
from helper_funcs.helpers_test import file_to_string


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
        self.assertTrue(
            datetime.now(timezone.utc) - alignment[0].created < timedelta(minutes=1),
            'alignment was not created within the last minute'
        )

    def test_save_alignment_adds_slug_value(self):
        """
        Tests that when an alignment is saved to the database it contains a slug value
        """
        Alignment.objects.create_alignment(self.name, self.data)
        alignment = Alignment.objects.all()
        self.assertTrue(len(alignment[0].slug) == 16, alignment[0].slug)
        slug_pattern = re.compile('^([a-zA-Z]|\d){16}$')
        self.assertTrue(re.match(slug_pattern, alignment[0].slug), alignment[0].slug)

    def test_saved_alignment_can_be_retrieved_with_slug(self):
        """
        Tests that when an alignment that is saved to the database can be retrieved by its slug value
        """
        alignment = Alignment.objects.create_alignment(self.name, self.data)
        slug = alignment.slug
        id = alignment.pk
        alignment_got = Alignment.objects.get(slug=slug)
        alignment_id = Alignment.objects.get(slug=slug).pk
        self.assertEqual('A. tha. SPA family alignment', alignment_got.name, alignment_got.name)
        self.assertEqual(id, alignment_id, alignment_id)

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
