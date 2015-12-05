from django.test import TestCase
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
        pass

    def test_form_validdation_for_invalid_characters(self):
        pass

    def test_form_validation_for_invalid_alignment(self):
        pass
