from django.test import TestCase
from django.utils.html import escape
from with_asserts.mixin import AssertHTMLMixin
import os
import io
from formalign.settings import BASE_DIR
from base.forms import QueryForm, EMPTY_ALIGNMENT_SUBMISSION_ERROR
from base.forms import parse_fasta

__author__ = 'Stefan Dieterle'


class IndexViewTestCase(TestCase):
    def test_index_view_basic(self):
        """
        Tests that index view returns a 200 response and uses the correct template
        :return:
        """
        response = self.client.get('/')
        self.assertTemplateUsed('base/index.html')
        self.assertEqual(response.status_code, 200)


class SeqDisplayTestCase(TestCase, AssertHTMLMixin):
    def test_display_page_uses_display_seq_template(self):
        """
        Tests that seq_display view returns a 200 response on a POST request and uses the correct template
        :return:
        """
        with open(os.path.join(BASE_DIR, 'test_data/short.fasta'), 'r') as alignment_file:
            input_seqs = alignment_file.read()
            response = self.client.post('/query-sequences/', {'align_input': input_seqs})
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed('base/query_display.html')

    def test_display_page_displays_query_seq(self):
        """
        Tests that seq_display displays the query on a valid POST request
        :return:
        """
        with open(os.path.join(BASE_DIR, 'test_data/short.fasta'), 'r') as alignment_file:
            input_seqs = alignment_file.read()
        response = self.client.post('/query-sequences/', {'align_input': input_seqs})
        with self.assertHTML(response, 'li[class=query_seq_meta]') as elems:
            self.assertEqual(elems[0].text,
                             'Short sequence1:',
                             'meta1: ' + format(elems[0].text)
                             )
            self.assertEqual(elems[1].text,
                             'Short sequence2:',
                             'meta2: ' + format(elems[1].text)
                             )
        with self.assertHTML(response, 'p[class=query_seq_display]') as elems:
            self.assertEqual(elems[0].text,
                             'MKERBGWAQ--QGKKPWRF--EEW',
                             'seq1: ' + elems[0].text
                             )
            self.assertNotIn(' ', elems[0].text)
            self.assertNotIn('\n', elems[0].text)
            self.assertEqual(elems[1].text,
                             'MKERBGWA-SYQGKKPWRFAQ-EW',
                             'seq2: ' + format(elems[1].text)
                             )

    def test_for_empty_input_renders_index_template(self):
        """
        Tests that submitting an empty alignment input renders the home page
        :return:
        """
        response = self.client.post('/query-sequences/', {'align_input': ''})
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed('base/index.html')

    def test_empty_input_errors_are_shown_on_home_page(self):
        """
        Tests that the correct error message is displayed on the home page on submission of an empty alignment
        :return:
        """
        response = self.client.post('/query-sequences/', {'align_input': ''})
        self.assertContains(response, escape(EMPTY_ALIGNMENT_SUBMISSION_ERROR))

    def test_for_empty_input_passes_form_to_template(self):
        """
        Tests that the form context is passed in response to an empty alignment submission
        :return:
        """
        response = self.client.post('/query-sequences/', {'align_input': ''})
        self.assertIsInstance(response.context['form'], QueryForm)

    def test_for_invalid_fasta_input_renders_index_template(self):
        self.fail('Incomplete Test')

    def test_invalid_fasta_input_errors_are_shown_on_home_page(self):
        self.fail('Incomplete Test')

    def test_for_invalid_fasta_input_passes_form_to_template(self):
        self.fail('Incomplete Test')

    def test_parse_fasta(self):
        """
        Tests that the parse_fasta function returns expected values with a valid fasta alignment
        :return:
        """
        with open(os.path.join(BASE_DIR, 'test_data/short.fasta'), 'r') as alignment_file:
            input_seqs = alignment_file.read()
        parsed = parse_fasta(io.StringIO(input_seqs))
        self.assertEqual(parsed[0]['meta'], 'Short sequence1')
        self.assertEqual(parsed[0]['seq'], 'MKERBGWAQ--QGKKPWRF--EEW')
        self.assertEqual(parsed[1]['meta'], 'Short sequence2')
        self.assertEqual(parsed[1]['seq'], 'MKERBGWA-SYQGKKPWRFAQ-EW')
