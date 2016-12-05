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

from base.forms import QueryForm
from base.forms import EMPTY_ERROR, FORMAT_ERROR, CHARACTER_ERROR, ALIGNMENT_ERROR, LESS_THAN_TWO_SEQS_ERROR
from helper_funcs.helpers_test import file_to_string
from Bio.Alphabet.IUPAC import ExtendedIUPACProtein, ExtendedIUPACDNA
from Bio.Alphabet import Gapped

__author__ = 'Stefan Dieterle'


class QueryFormTest(TestCase):

    def validation(self, error_text, input_file='', seq_type='protein'):
        """
        Performs validation test for invalid forms, takes a user alignment
        input, asserts form.is_valid as false and checks the error
        :param:
        :return:
        """
        if input_file:
            input_seqs = file_to_string(input_file)
        else:
            input_seqs = ''
        form = QueryForm(data={'align_input': input_seqs, 'seq_type': seq_type})
        self.assertFalse(form.is_valid())
        self.assertEqual(
            form.errors['align_input'],
            [error_text], format(form.errors['align_input'])
        )

    def test_form_renders_seq_text_input_widget(self):
        """
        Tests correct rendering of form input text area
        :return:
        """
        form = QueryForm()
        field = form['align_input']
        self.assertIn('Paste in your alignment:<br>(FASTA, clustalw, stockholm or phylip)', field.label_tag())
        self.assertIn('placeholder="Alignment (FASTA, clustalw, stockholm or phylip)"', form.as_p())
        self.assertIn('class="form-control"', form.as_p())

    def test_form_renders_sequence_type_radio_buttons(self):
        """
        Tests correct rendering of form sequence type radio buttons
        :return:
        """
        form = QueryForm()
        field = form['seq_type']
        self.assertIn('Input sequence type:', field.label_tag())
        self.assertIn('<input id="id_seq_type_0" name="seq_type" type="radio" value="protein"', form.as_p())
        self.assertIn(
            '<input checked="checked" id="id_seq_type_1" name="seq_type" type="radio" value="DNA"', form.as_p()
        )

    def test_form_renders_consensus_type_radio_buttons(self):
        """
        Tests correct rendering of form consensus type radio buttons
        :return:
        """
        form = QueryForm()
        field = form['cons_type']
        self.assertIn('Consensus type:', field.label_tag())
        self.assertIn('<input id="id_cons_type_1" name="cons_type" type="radio" value="subs_matrix"', form.as_p())
        self.assertIn(
            '<input checked="checked" id="id_cons_type_0" name="cons_type" type="radio" value="identity"', form.as_p()
        )

    def test_form_validation_for_blank_items(self):
        """
        Tests error on submission of empty form
        :return:
        """
        self.validation(EMPTY_ERROR)
        self.validation(EMPTY_ERROR, '', 'DNA')

    def test_form_validation_for_invalid_fasta(self):
        """
        Tests error on invalid FASTA (no '>' as first character)
        :return:
        """
        self.validation(FORMAT_ERROR, 'protein_invalid_fasta.fasta')
        self.validation(FORMAT_ERROR, 'DNA_invalid_fasta.fasta', 'DNA')

    def test_form_validation_for_invalid_characters(self):
        """
        Tests error on invalid characters in FASTA sequence for protein sequences
        :return:
        """
        self.validation(CHARACTER_ERROR + 'sequence1', 'protein_invalid_characters.fasta')
        self.validation(CHARACTER_ERROR + 'sequence1', 'DNA_invalid_characters.fasta', 'DNA')

    def test_form_validation_for_invalid_alignment(self):
        """
        Tests error on invalid alignment when sequences have differing lengths
        :return:
        """
        self.validation(ALIGNMENT_ERROR, 'protein_invalid_alignment.fasta')
        self.validation(ALIGNMENT_ERROR, 'DNA_invalid_alignment.fasta', 'DNA')

    def test_form_validation_for_too_few_sequences_in_alignment(self):
        """
        Tests error on invalid alignment when sequences have differing lengths
        :return:
        """
        self.validation(LESS_THAN_TWO_SEQS_ERROR, 'protein_too_few_sequences.fasta')
        self.validation(LESS_THAN_TWO_SEQS_ERROR, 'DNA_too_few_sequences.fasta', 'DNA')

    def test_form_validation_returns_correct_seqrecord_alphabet(self):
        """
        Tests that seqrecords get correct alphabet according to user input of seq_type
        :return:
        """
        input_seqs = file_to_string('protein.fasta')
        form = QueryForm(data={'align_input': input_seqs, 'seq_type': 'protein', 'cons_type': 'identity'})
        self.assertTrue(form.is_valid())
        for f in form.cleaned_data['align_input']:
            self.assertEqual(
                "Gapped(ExtendedIUPACProtein(), '-')",
                str(f.seq.alphabet),
                Gapped(ExtendedIUPACProtein()).letters
            )

        input_seqs = file_to_string('DNA.fasta')
        form = QueryForm(data={'align_input': input_seqs, 'seq_type': 'DNA', 'cons_type': 'identity'})
        self.assertTrue(form.is_valid())
        for f in form.cleaned_data['align_input']:
            self.assertEqual(
                "Gapped(ExtendedIUPACDNA(), '-')",
                str(f.seq.alphabet),
                Gapped(ExtendedIUPACDNA()).letters
            )
