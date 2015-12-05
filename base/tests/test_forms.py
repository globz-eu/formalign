from django.test import TestCase
import os
from formalign.settings import BASE_DIR
from base.forms import QueryForm

__author__ = 'Stefan Dieterle'


class QueryFormTest(TestCase):

    def test_form_renders_seq_text_input(self):
        """
        Tests correct rendering of form elements
        :return:
        """
        form = QueryForm()
        field = form['align_input']
        self.assertIn('Paste in your alignment in FASTA format:', field.label_tag())
        self.assertIn('placeholder="FASTA alignment"', form.as_p())
        self.assertIn('class="form-control"', form.as_p())

    def test_form_validation_for_blank_items(self):
        """
        Tests error on submission of empty form
        :return:
        """
        form = QueryForm(data={'align_input': ''})
        self.assertFalse(form.is_valid())
        self.assertEqual(
            form.errors['align_input'],
            ['Please submit an alignment']
        )

    def test_form_validation_for_invalid_fasta(self):
        """
        Tests error on invalid FASTA (no '>' as first character)
        :return:
        """
        with open(os.path.join(BASE_DIR, 'test_data/short_invalid_fasta.fasta'), 'r') as alignment_file:
            input_seqs = alignment_file.read()
        form = QueryForm(data={'align_input': input_seqs})
        self.assertFalse(form.is_valid())
        self.assertEqual(
            form.errors['align_input'],
            ['Sequence is not FASTA compliant, no ">" as first character']
        )

    def test_form_validation_for_invalid_characters(self):
        """
        Tests error on invalid characters in FASTA sequence for protein sequences
        :return:
        """
        with open(os.path.join(BASE_DIR, 'test_data/short_invalid_characters.fasta'), 'r') as alignment_file:
            input_seqs = alignment_file.read()
        form = QueryForm(data={'align_input': input_seqs})
        self.assertFalse(form.is_valid())
        self.assertEqual(
            form.errors['align_input'],
            ['Invalid character in sequence: Short sequence3']
        )

    def test_form_validation_for_invalid_alignment(self):
        """
        Tests error on invalid alignment when sequences have differing lengths
        :return:
        """
        with open(os.path.join(BASE_DIR, 'test_data/short_invalid_alignment.fasta'), 'r') as alignment_file:
            input_seqs = alignment_file.read()
        form = QueryForm(data={'align_input': input_seqs})
        self.assertFalse(form.is_valid())
        self.assertEqual(
            form.errors['align_input'],
            ['Alignment invalid, sequences have different lengths']
        )
