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
from datetime import datetime, timedelta, timezone

from Bio.Alphabet import Gapped
from Bio.Alphabet.IUPAC import ExtendedIUPACProtein
from base.models import Alignment, Seqrecord
from base.tasks import clean_alignments
from django.test import TestCase
from django.test.utils import override_settings
from formalign.settings import CLEAN_OLDER
from helper_funcs.bio.helpers import parse_fasta_alignment
from helper_funcs.helpers_test import file_to_string

__author__ = 'Stefan Dieterle'


class CleanAlignmentsTestCase(TestCase):
    """
    Tests the clean_alignments task
    """

    def setUp(self):
        self.name = 'A. tha. SPA family alignment'
        align_input = io.StringIO(file_to_string('spa_protein_alignment.fasta'))
        self.data = parse_fasta_alignment(align_input)
        alphabet = Gapped(ExtendedIUPACProtein())
        for a in self.data:
            a.seq.alphabet = alphabet
        self.data._alphabet = alphabet
        self.older = int(CLEAN_OLDER)

    @override_settings(CELERY_EAGER_PROPAGATES_EXCEPTIONS=True, CELERY_ALWAYS_EAGER=True)
    def test_clean_alignments_removes_old_alignments(self):
        """
        Tests that an alignment older than 7 days is removed from the database when clean_alignments is run
        """
        alignment = Alignment.objects.create_alignment(self.name, self.data)
        alignment.created = datetime.now(timezone.utc) - timedelta(days=self.older)
        alignment.save()
        modified = Alignment.objects.get(pk=alignment.pk)
        self.assertTrue(datetime.now(timezone.utc) - modified.created > timedelta(days=self.older), modified.created)
        response = self.client.get('/align-display/%s' % alignment.slug)
        self.assertEqual(200, response.status_code, response.status_code)
        result = clean_alignments.delay()
        r = result.get()
        self.assertEqual(['A. tha. SPA family alignment'], r, r)
        response = self.client.get('/align-display/%s' % alignment.slug)
        self.assertEqual(404, response.status_code, response.status_code)
        try:
            old = Alignment.objects.get(pk=alignment.pk)
            self.assertFalse(old.name == self.name, '%s was not removed' % old.name)
        except Alignment.DoesNotExist as error:
            self.assertEqual(str(error), 'Alignment matching query does not exist.', error)

    @override_settings(CELERY_EAGER_PROPAGATES_EXCEPTIONS=True, CELERY_ALWAYS_EAGER=True)
    def test_clean_alignments_removes_sequences_associated_to_old_alignments(self):
        """
        Tests that all sequences associated to an alignment older than 7 days are removed from the database when
        clean_alignments is run
        """
        alignment = Alignment.objects.create_alignment(self.name, self.data)
        alignment.created = datetime.now(timezone.utc) - timedelta(days=self.older)
        alignment.save()
        modified = Alignment.objects.get(pk=alignment.pk)
        self.assertTrue(datetime.now(timezone.utc) - modified.created > timedelta(days=self.older), modified.created)
        seq_ids = [m.pk for m in modified.seqs.all()]
        result = clean_alignments.delay()
        r = result.get()
        self.assertEqual(['A. tha. SPA family alignment'], r, r)
        for s in seq_ids:
            try:
                old_seq = Seqrecord.objects.get(pk=s)
                self.assertFalse(old_seq.pk in seq_ids, 'sequence %s was not removed' % old_seq.pk)
            except Seqrecord.DoesNotExist as error:
                self.assertEqual(str(error), 'Seqrecord matching query does not exist.', error)

    @override_settings(CELERY_EAGER_PROPAGATES_EXCEPTIONS=True, CELERY_ALWAYS_EAGER=True)
    def test_clean_alignments_leaves_recent_alignments_untouched(self):
        """
        Tests that an alignment not older than 7 days is not removed from the database when clean_alignments is run
        """
        alignment = Alignment.objects.create_alignment(self.name, self.data)
        alignment.created = datetime.now(timezone.utc) - timedelta(days=self.older) + timedelta(hours=1)
        alignment.save()
        modified = Alignment.objects.get(pk=alignment.pk)
        self.assertTrue(datetime.now(timezone.utc) - modified.created < timedelta(days=self.older), modified.created)
        result = clean_alignments.delay()
        r = result.get()
        self.assertEqual([], r, r)
        try:
            new = Alignment.objects.get(name=self.name)
            self.assertTrue(new.name == self.name, '%s was removed' % new.name)
        except Alignment.DoesNotExist as error:
            self.assertFalse(str(error) == 'Alignment matching query does not exist.', error)

    @override_settings(CELERY_EAGER_PROPAGATES_EXCEPTIONS=True, CELERY_ALWAYS_EAGER=True)
    def test_clean_alignments_leaves_sequences_associated_to_recent_alignments_untouched(self):
        """
        Tests that all sequences associated to an alignment not older than 7 days are not removed from the database when
        clean_alignments is run
        """
        alignment = Alignment.objects.create_alignment(self.name, self.data)
        alignment.created = datetime.now(timezone.utc) - timedelta(days=self.older) + timedelta(hours=1)
        alignment.save()
        modified = Alignment.objects.get(pk=alignment.pk)
        self.assertTrue(datetime.now(timezone.utc) - modified.created < timedelta(days=self.older), modified.created)
        seq_ids = [m.pk for m in modified.seqs.all()]
        result = clean_alignments.delay()
        r = result.get()
        self.assertEqual([], r, r)
        for s in seq_ids:
            try:
                new_seq = Seqrecord.objects.get(pk=s)
                self.assertTrue(new_seq.pk in seq_ids, 'sequence %s was removed' % new_seq.pk)
            except Seqrecord.DoesNotExist as error:
                self.assertFalse(str(error) == 'Seqrecord matching query does not exist.', error)

    @override_settings(CELERY_EAGER_PROPAGATES_EXCEPTIONS=True, CELERY_ALWAYS_EAGER=True)
    def test_clean_alignments_cleans_only_old_alignemnts(self):
        """
        Tests that only alignments older than 7 days are removed from the database when clean_alignments is run
        """
        alignment_old = Alignment.objects.create_alignment(self.name, self.data)
        alignment_new = Alignment.objects.create_alignment(self.name, self.data)
        alignment_old.created = datetime.now(timezone.utc) - timedelta(days=self.older)
        alignment_old.save()
        modified_old = Alignment.objects.get(pk=alignment_old.pk)
        self.assertTrue(
            datetime.now(timezone.utc) - modified_old.created > timedelta(days=self.older),
            modified_old.created
        )

        result = clean_alignments.delay()
        r = result.get()
        self.assertEqual(['A. tha. SPA family alignment'], r, r)

        try:
            new = Alignment.objects.get(pk=alignment_new.pk)
            self.assertTrue(new.name == self.name, '%s was not removed' % new.name)
        except Alignment.DoesNotExist as error:
            self.assertFalse(str(error) == 'Alignment matching query does not exist.', error)
        try:
            old = Alignment.objects.get(pk=alignment_old.pk)
            self.assertFalse(old.name == self.name, '%s was removed' % old.name)
        except Alignment.DoesNotExist as error:
            self.assertEqual(str(error), 'Alignment matching query does not exist.', error)

    @override_settings(CELERY_EAGER_PROPAGATES_EXCEPTIONS=True, CELERY_ALWAYS_EAGER=True)
    def test_clean_alignments_cleans_only_old_alignemnts_sequences(self):
        """
        Tests that only alignments older than 7 days are removed from the database when clean_alignments is run
        """
        seq_ids_old = []
        for i in range(2):
            alignment_old = Alignment.objects.create_alignment(self.name, self.data)
            alignment_old.created = datetime.now(timezone.utc) - timedelta(days=self.older)
            alignment_old.save()
            modified_old = Alignment.objects.get(pk=alignment_old.pk)
            seq_ids_old.extend([m.pk for m in modified_old.seqs.all()])
        alignment_new = Alignment.objects.create_alignment(self.name, self.data)
        seq_ids_new = [a.pk for a in alignment_new.seqs.all()]

        result = clean_alignments.delay()
        r = result.get()
        self.assertEqual(['A. tha. SPA family alignment', 'A. tha. SPA family alignment'], r, r)

        for s in seq_ids_new:
            try:
                new_seq = Seqrecord.objects.get(pk=s)
                self.assertTrue(new_seq.pk in seq_ids_new, 'sequence %s was removed' % new_seq.pk)
            except Seqrecord.DoesNotExist as error:
                self.assertFalse(str(error) == 'Seqrecord matching query does not exist.', error)
        for s in seq_ids_old:
            try:
                old_seq = Seqrecord.objects.get(pk=s)
                self.assertFalse(old_seq.pk in seq_ids_old, 'sequence %s was not removed' % old_seq.pk)
            except Seqrecord.DoesNotExist as error:
                self.assertEqual(str(error), 'Seqrecord matching query does not exist.', error)
