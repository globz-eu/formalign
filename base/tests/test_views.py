import io
import re

from django.test import TestCase
from django.utils.html import escape
from with_asserts.mixin import AssertHTMLMixin

from base.forms import QueryForm
from base.forms import EMPTY_ERROR, FORMAT_ERROR, CHARACTER_ERROR, ALIGNMENT_ERROR, LESS_THAN_TWO_SEQS_ERROR
from helper_funcs.helpers_bio import parse_fasta_alignment
from helper_funcs.helpers_test import file_to_string

from Bio.Alphabet.IUPAC import ExtendedIUPACProtein, ExtendedIUPACDNA
from Bio.Alphabet import Gapped

from base.models import save_alignment_to_db, get_multipleseqalignment_object_from_db

__author__ = 'Stefan Dieterle'


class IndexViewTestCase(TestCase, AssertHTMLMixin):
    """
    Tests for index view
    """

    def test_index_view_basic(self):
        """
        Tests that index view returns a 200 response and uses the correct template
        :return:
        """
        response = self.client.get('/')
        self.assertTemplateUsed(response, 'base/index.html')
        self.assertEqual(response.status_code, 200)

    def test_redirect_to_seqdisplay_on_post(self):
        """
        Tests that valid POST request on index page redirects to /query-sequences/
        :return:
        """
        input_seqs = file_to_string('spa_protein_alignment.fasta')
        response = self.client.post('/', {'align_input': input_seqs, 'seq_type': 'Protein'})
        self.assertTrue('/query-sequences/' in response.url)

    def test_alignment_is_saved_on_post(self):
        """
        Tests that alignment is saved on valid POST request to index
        :return:
        """
        input_seqs = file_to_string('spa_protein_alignment.fasta')
        response = self.client.post('/', {'align_input': input_seqs, 'seq_type': 'Protein'})
        print(response.url)
        pk = re.match(r'^/query-sequences/(?P<align_id>\d+)/', response.url).group('align_id')
        alignment = get_multipleseqalignment_object_from_db(pk)
        # self.assertEqual('bla', alignment.seqs, alignment.seqs)
        self.assertEqual(
                [seq.id for seq in alignment],
                ['NP_175717', 'NP_683567', 'NP_182157', 'NP_192849']
        )


class SeqDisplayTestCase(TestCase, AssertHTMLMixin):
    """
    Tests for sequence display
    """

    def setUp(self):
        name = 'A. tha. SPA family protein alignment'
        align_input = io.StringIO(file_to_string('spa_protein_alignment.fasta'))
        data = parse_fasta_alignment(align_input)
        for d in data:
            d.seq.alphabet = Gapped(ExtendedIUPACProtein())
        align_pk = save_alignment_to_db(name, data)
        self.response_prot = self.client.get('/query-sequences/' + str(align_pk) + '/')
        name = 'A. tha. SPA family DNA alignment'
        align_input = io.StringIO(file_to_string('spa_cds_alignment.fasta'))
        data = parse_fasta_alignment(align_input)
        for d in data:
            d.seq.alphabet = Gapped(ExtendedIUPACDNA())
        align_pk = save_alignment_to_db(name, data)
        self.response_dna = self.client.get('/query-sequences/' + str(align_pk) + '/')

    def test_display_page_uses_display_seq_template(self):
        """
        Tests that seq_display view returns a 200 response on a POST request and uses the correct template
        :return:
        """
        input_seqs = file_to_string('protein.fasta')
        response = self.client.post('/', {'align_input': input_seqs, 'seq_type': 'DNA'})
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, 'base/index.html')

    def test_display_page_displays_sequence_type(self):
        """
        Tests that seq_display displays the selected sequence type on a valid POST request
        :return:
        """
        response = self.response_prot
        with self.assertHTML(response, 'h2[class="query_seq_type color-complement-4"]') as elem:
            self.assertEqual(elem[0].text,
                             'Protein sequences:',
                             format(elem[0].text)
                             )

        response = self.response_dna
        with self.assertHTML(response, 'h2[class="query_seq_type color-complement-4"]') as elem:
            self.assertEqual(elem[0].text,
                             'DNA sequences:',
                             format(elem[0].text)
                             )

    def test_display_page_displays_protein_query_seq(self):
        """
        Tests that seq_display displays the query on a valid POST request
        :return:
        """
        response = self.response_prot
        with self.assertHTML(response, 'h3[class="query_seq_meta bg-color-body"]') as elems:
            self.assertEqual('NP_175717 NP_175717.1 SPA1-related 4 protein [Arabidopsis thaliana].:',
                             elems[0].text,
                             'meta1: ' + format(elems[0].text)
                             )
            self.assertEqual('NP_683567 NP_683567.1 protein SPA1-related 3 [Arabidopsis thaliana].:',
                             elems[1].text,
                             'meta2: ' + format(elems[1].text)
                             )
        with self.assertHTML(response, 'p[class=query_seq_display]') as elems:
            self.assertEqual(elems[0].text,
                             '-' * 80,
                             'seq1: ' + elems[0].text
                             )
            self.assertNotIn(' ', elems[0].text)
            self.assertNotIn('\n', elems[0].text)
            self.assertEqual(elems[1].text,
                             '-' * 80,
                             'seq2: ' + format(elems[1].text)
                             )

    def test_display_page_displays_sequence_with_less_than_80_residues_per_line(self):
        response = self.response_prot
        with self.assertHTML(response, 'p[class=query_seq_display]') as elems:
            for elem in elems:
                self.assertTrue(len(elem.text) <= 80)

    def test_display_page_displays_consensus(self):
        """
        Tests that seq_display displays the consensus sequence on a valid POST request
        :return:
        """
        response = self.response_prot
        with self.assertHTML(response, 'h3[class="query_seq_meta bg-color-body"]') as elems:
            self.assertEqual(elems[-1].text,
                             'Consensus:',
                             'consensus meta: ' + format(elems[0].text)
                             )
        cons_seq = file_to_string('consensus.txt')
        with self.assertHTML(response, 'p[class="query_seq_display consensus_seq_display"]') as elems:
            self.assertEqual(elems[0].text,
                             cons_seq[:80],
                             'consensus seq: ' + elems[0].text
                             )
            self.assertNotIn(' ', elems[0].text)
            self.assertNotIn('\n', elems[0].text)

    def test_parse_fasta_alignment(self):
        """
        Tests that the parse_fasta function returns expected values with a valid fasta alignment
        :return:
        """
        input_seqs = file_to_string('protein.fasta')
        parsed = parse_fasta_alignment(io.StringIO(input_seqs))
        self.assertEqual(parsed[0].description, 'sequence1')
        self.assertEqual(parsed[0].seq, 'MKERBGWAQ--QGKKPWRF--EEW')
        self.assertEqual(parsed[1].description, 'sequence2')
        self.assertEqual(parsed[1].seq, 'MKERBGWA-SYQGKKPWRFAQ-EW')

    def test_display_page_uses_display_seq_template_on_GET(self):
        """
        Tests that seq_display view returns a 200 response on a GET request and uses the correct template
        :return:
        """
        name = 'A. tha. SPA family alignment'
        align_input = io.StringIO(file_to_string('spa_protein_alignment.fasta'))
        data = parse_fasta_alignment(align_input)
        for d in data:
            d.seq.alphabet = Gapped(ExtendedIUPACProtein())
        save = save_alignment_to_db(name, data)
        response = self.client.get('/query-sequences/' + str(save) + '/')
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, 'base/query_display.html')


class SeqDisplayInvalidInput(TestCase):
    """
    Tests for invalid alignment submission
    """

    def response_for_invalid_post_request(self, input_file='', seq_type='Protein'):
        """
        Creates a response from a POST request to /query-sequences/ with an invalid alignment
        :param input_file: file containing invalid alignment
        :return: response
        """
        if input_file:
            input_seqs = file_to_string(input_file)
        else:
            input_seqs = ''
        response = self.client.post('/', {'align_input': input_seqs, 'seq_type': seq_type})
        return response

    def invalid_input_renders_index_template(self, input_file='', seq_type='Protein'):
        """
        Tests that submitting an invalid alignment input renders the home page
        :param input_file: file containing invalid alignment
        :return:
        """
        response = self.response_for_invalid_post_request(input_file, seq_type)
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, 'base/index.html')

    def invalid_input_errors_are_shown_on_home_page(self, error_text, input_file='', seq_type='Protein'):
        """
        Tests that the correct error message is displayed on the home page on submission of an invalid alignment
        :param error_text: expected error message
        :param input_file: file containing invalid alignment
        :return:
        """
        response = self.response_for_invalid_post_request(input_file, seq_type)
        self.assertContains(response, escape(error_text))

    def invalid_input_passes_form_to_template(self, input_file='', seq_type='Protein'):
        """
        Tests that the form context is passed in response to an invalid alignment submission
        :param input_file: file containing invalid alignment
        :return:
        """
        response = self.response_for_invalid_post_request(input_file, seq_type)
        self.assertIsInstance(response.context['form'], QueryForm)

    def test_for_empty_input_renders_index_template(self):
        """
        Tests that submitting an empty alignment input renders the home page
        :return:
        """
        self.invalid_input_renders_index_template('')
        self.invalid_input_renders_index_template('', 'DNA')

    def test_empty_input_errors_are_shown_on_home_page(self):
        """
        Tests that the correct error message is displayed on the home page on submission of an empty alignment
        :return:
        """
        self.invalid_input_errors_are_shown_on_home_page(EMPTY_ERROR)
        self.invalid_input_errors_are_shown_on_home_page(EMPTY_ERROR, '', 'DNA')

    def test_for_empty_input_passes_form_to_template(self):
        """
        Tests that the form context is passed in response to an empty alignment submission
        :return:
        """
        self.invalid_input_passes_form_to_template()
        self.invalid_input_passes_form_to_template(seq_type='DNA')

    def test_for_invalid_fasta_input_renders_index_template(self):
        """
        Tests that submitting an empty alignment input renders the home page
        :return:
        """
        self.invalid_input_renders_index_template('protein_invalid_fasta.fasta')
        self.invalid_input_renders_index_template('DNA_invalid_fasta.fasta', 'DNA')

    def test_invalid_fasta_input_errors_are_shown_on_home_page(self):
        """
        Tests that the correct error message is displayed on the home page on submission of invalid FASTA
        :return:
        """
        self.invalid_input_errors_are_shown_on_home_page(FORMAT_ERROR, 'protein_invalid_fasta.fasta')
        self.invalid_input_errors_are_shown_on_home_page(FORMAT_ERROR, 'DNA_invalid_fasta.fasta', 'DNA')

    def test_for_invalid_fasta_input_passes_form_to_template(self):
        """
        Tests that the form context is passed in response to an invalid FASTA submission
        :return:
        """
        self.invalid_input_passes_form_to_template('protein_invalid_fasta.fasta')
        self.invalid_input_passes_form_to_template('DNA_invalid_fasta.fasta', 'DNA')

    def test_for_invalid_character_input_renders_index_template(self):
        """
        Tests that submitting an empty alignment input renders the home page
        :return:
        """
        self.invalid_input_renders_index_template('protein_invalid_characters.fasta')
        self.invalid_input_renders_index_template('DNA_invalid_fasta.fasta', 'DNA')

    def test_invalid_character_input_errors_are_shown_on_home_page(self):
        """
        Tests that the correct error message is displayed on the home page on submission of invalid FASTA
        :return:
        """
        self.invalid_input_errors_are_shown_on_home_page(CHARACTER_ERROR, 'protein_invalid_characters.fasta')
        self.invalid_input_errors_are_shown_on_home_page(CHARACTER_ERROR, 'DNA_invalid_characters.fasta', 'DNA')

    def test_for_invalid_character_input_passes_form_to_template(self):
        """
        Tests that the form context is passed in response to an invalid FASTA submission
        :return:
        """
        self.invalid_input_passes_form_to_template('protein_invalid_characters.fasta')
        self.invalid_input_passes_form_to_template('DNA_invalid_characters.fasta', 'DNA')

    def test_for_invalid_alignment_input_renders_index_template(self):
        """
        Tests that submitting an empty alignment input renders the home page
        :return:
        """
        self.invalid_input_renders_index_template('protein_invalid_alignment.fasta')
        self.invalid_input_renders_index_template('DNA_invalid_alignment.fasta', 'DNA')

    def test_invalid_alignment_input_errors_are_shown_on_home_page(self):
        """
        Tests that the correct error message is displayed on the home page on submission of invalid FASTA
        :return:
        """
        self.invalid_input_errors_are_shown_on_home_page(ALIGNMENT_ERROR, 'protein_invalid_alignment.fasta')
        self.invalid_input_errors_are_shown_on_home_page(ALIGNMENT_ERROR, 'DNA_invalid_alignment.fasta', 'DNA')

    def test_for_invalid_alignment_input_passes_form_to_template(self):
        """
        Tests that the form context is passed in response to an invalid FASTA submission
        :return:
        """
        self.invalid_input_passes_form_to_template('protein_invalid_alignment.fasta')
        self.invalid_input_passes_form_to_template('DNA_invalid_alignment.fasta', 'DNA')

    def test_for_too_few_sequences_input_renders_index_template(self):
        """
        Tests that submitting an alignment with only one sequence input renders the home page
        :return:
        """
        self.invalid_input_renders_index_template('protein_too_few_sequences.fasta')
        self.invalid_input_renders_index_template('DNA_too_few_sequences.fasta', 'DNA')

    def test_too_few_sequences_input_errors_are_shown_on_home_page(self):
        """
        Tests that the correct error message is displayed on the home page on submission of invalid FASTA
        :return:
        """
        self.invalid_input_errors_are_shown_on_home_page(LESS_THAN_TWO_SEQS_ERROR, 'protein_too_few_sequences.fasta')
        self.invalid_input_errors_are_shown_on_home_page(LESS_THAN_TWO_SEQS_ERROR, 'DNA_too_few_sequences.fasta', 'DNA')

    def test_for_too_few_sequences_input_passes_form_to_template(self):
        """
        Tests that the form context is passed in response to an invalid FASTA submission
        :return:
        """
        self.invalid_input_passes_form_to_template('protein_too_few_sequences.fasta')
        self.invalid_input_passes_form_to_template('DNA_too_few_sequences.fasta', 'DNA')


class AlignDisplayTestCase(TestCase, AssertHTMLMixin):
    """
    Tests for alignment display
    """
    def setUp(self):
        """
        Creates a response from a GET request to /align-display/ with an alignment pk
        :param input_file: file containing alignment
        :return: response
        """
        name = 'A. tha. SPA family protein alignment'
        align_input = io.StringIO(file_to_string('spa_protein_alignment.fasta'))
        data = parse_fasta_alignment(align_input)
        for d in data:
            d.seq.alphabet = Gapped(ExtendedIUPACProtein())
        align_pk = save_alignment_to_db(name, data)
        self.response = self.client.get('/align-display/' + str(align_pk) + '/')

    def test_align_display_page_uses_align_display_seq_template(self):
        """
        Tests that align_display view returns a 200 response on a GET request and uses the correct template
        :return:
        """
        self.assertEqual(self.response.status_code, 200)
        self.assertTemplateUsed(self.response, 'base/align_display.html')

    def test_align_display_page_displays_default_formatted_protein_alignment(self):
        """
        Tests that align_display displays an alignment in the default format on a valid GET request (line of 80
        characters in blocks of 10 characters)
        :return:
        """
        with self.assertHTML(self.response, 'div[class="al_ln"]') as elems:
            self.assertEqual(len([e for e in elems[0].findall('div') if e.attrib['class'] != 'sep']),
                             81, elems[0].getchildren()[0].text
                             )
            self.assertEqual(len([e for e in elems[0].findall('div') if e.attrib['class'] == 'sep']),
                             8, elems[0].getchildren()[0].text
                             )

    def test_align_display_page_displays_correct_protein_alignment_sequence(self):
        """
        Tests that align_display displays an alignment with correct sequences
        :return:
        """
        expected_seqs = file_to_string('spa_protein_alignment.fasta')
        align_expected = io.StringIO(expected_seqs)
        alignment = parse_fasta_alignment(align_expected)
        seq = alignment[3].seq
        with self.assertHTML(self.response, 'div[class="al_ln"]') as elems:
            for i, e in enumerate(
                    [
                        elem for elem in elems[3].findall('div')
                        if elem.attrib['class'] not in ['sep', 'al_el seq_id']
                        ]
            ):
                self.assertEqual(e.text, seq[i], 'e.text: ' + format(e.attrib))
