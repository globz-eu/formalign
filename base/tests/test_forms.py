from django.test import TestCase
import os
from formalign.settings import BASE_DIR
from base.forms import QueryForm
from base.forms import EMPTY_ERROR, FASTA_ERROR, CHARACTER_ERROR, ALIGNMENT_ERROR, LESS_THAN_TWO_SEQS_ERROR
from helper_funcs.helpers_test import file_to_string

__author__ = 'Stefan Dieterle'


class QueryFormTest(TestCase):

    def validation(self, error_text, input_file=''):
        """
        Performs validation test for invalid forms, takes a user alignment input, asserts form.is_valid as false and
        checks the error
        :param:
        :return:
        """
        if input_file:
            input_seqs = file_to_string(input_file)
        else:
            input_seqs = ''
        form = QueryForm(data={'align_input': input_seqs})
        self.assertFalse(form.is_valid())
        self.assertEqual(
            form.errors['align_input'],
            [error_text]
        )

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
        self.validation(EMPTY_ERROR)

    def test_form_validation_for_invalid_fasta(self):
        """
        Tests error on invalid FASTA (no '>' as first character)
        :return:
        """
        self.validation(FASTA_ERROR, 'short_invalid_fasta.fasta')

    def test_form_validation_for_invalid_characters(self):
        """
        Tests error on invalid characters in FASTA sequence for protein sequences
        :return:
        """
        self.validation(CHARACTER_ERROR + 'Short sequence3', 'short_invalid_characters.fasta')

    def test_form_validation_for_invalid_alignment(self):
        """
        Tests error on invalid alignment when sequences have differing lengths
        :return:
        """
        self.validation(ALIGNMENT_ERROR, 'short_invalid_alignment.fasta')

    def test_form_validation_for_too_few_sequences_in_alignment(self):
        """
        Tests error on invalid alignment when sequences have differing lengths
        :return:
        """
        self.validation(LESS_THAN_TWO_SEQS_ERROR, 'short_too_few_sequences.fasta')
