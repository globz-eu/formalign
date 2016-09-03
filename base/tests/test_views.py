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
import re
import time
import pickle
import os

from django.test import TestCase
from django.utils.html import escape
from with_asserts.mixin import AssertHTMLMixin

from base.forms import QueryForm
from base.forms import EMPTY_ERROR, FORMAT_ERROR, CHARACTER_ERROR, ALIGNMENT_ERROR, LESS_THAN_TWO_SEQS_ERROR
from helper_funcs.helpers_bio import parse_fasta_alignment, consensus_add
from helper_funcs.helpers_test import file_to_string
from helper_funcs.helpers_format import split_lines, split_lines_in_blocks, annotate

from Bio.Alphabet.IUPAC import ExtendedIUPACProtein, ExtendedIUPACDNA
from Bio.Alphabet import Gapped
from Bio.Align import AlignInfo

from base.models import Alignment

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

        # Does not work with Jinja2
        # self.assertTemplateUsed(response, 'base/index.html')

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
        pk = re.match(r'^/query-sequences/(?P<align_id>\d+)/', response.url).group('align_id')
        alignment = Alignment.objects.get_alignment(pk)
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
        align = Alignment.objects.create_alignment(name, data)
        self.response_prot = self.client.get('/query-sequences/' + str(align.id) + '/')
        name = 'A. tha. SPA family DNA alignment'
        align_input = io.StringIO(file_to_string('spa_cds_alignment.fasta'))
        data = parse_fasta_alignment(align_input)
        for d in data:
            d.seq.alphabet = Gapped(ExtendedIUPACDNA())
        align = Alignment.objects.create_alignment(name, data)
        self.response_dna = self.client.get('/query-sequences/' + str(align.id) + '/')

    def test_display_page_uses_display_seq_template(self):
        """
        Tests that seq_display view returns a 200 response on a POST request and uses the correct template
        :return:
        """
        input_seqs = file_to_string('protein.fasta')
        response = self.client.post('/', {'align_input': input_seqs, 'seq_type': 'DNA'})
        self.assertEqual(response.status_code, 200)

        # Does not work with Jinja2
        # self.assertTemplateUsed(response, 'base/index.html')

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

    def test_display_page_displays_protein_sequence_with_less_than_80_residues_per_line(self):
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
                             'consensus 70%:',
                             'consensus meta: ' + format(elems[0].text)
                             )
        cons_seq = file_to_string('consensus.txt')
        with self.assertHTML(response, 'div[class="query_seq bg-color-body"]') as elems:
            self.assertEqual(elems[-1].findall('p')[0].text,
                             cons_seq[:80],
                             'consensus seq: ' + elems[-1].findall('p')[0].text
                             )
            self.assertNotIn(' ', elems[-1].findall('p')[0].text)
            self.assertNotIn('\n', elems[-1].findall('p')[0].text)

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
        save = Alignment.objects.create_alignment(name, data)
        response = self.client.get('/query-sequences/' + str(save.id) + '/')
        self.assertEqual(response.status_code, 200)

        # Does not work with Jinja2
        # self.assertTemplateUsed(response, 'base/query_display.html')


class SeqDisplayInvalidInput(TestCase):
    """
    Tests for invalid alignment submission
    """

    files = [
        {'file': '', 'seq_type': 'DNA', 'error_text': EMPTY_ERROR},
        {'file': '', 'seq_type': 'Protein', 'error_text': EMPTY_ERROR},
        {'file': 'DNA_invalid_fasta.fasta', 'seq_type': 'DNA', 'error_text': FORMAT_ERROR},
        {'file': 'protein_invalid_fasta.fasta', 'seq_type': 'Protein', 'error_text': FORMAT_ERROR},
        {'file': 'DNA_invalid_characters.fasta', 'seq_type': 'DNA', 'error_text': CHARACTER_ERROR},
        {'file': 'protein_invalid_characters.fasta', 'seq_type': 'Protein', 'error_text': CHARACTER_ERROR},
        {'file': 'DNA_invalid_alignment.fasta', 'seq_type': 'DNA', 'error_text': ALIGNMENT_ERROR},
        {'file': 'protein_invalid_alignment.fasta', 'seq_type': 'Protein', 'error_text': ALIGNMENT_ERROR},
        {'file': 'DNA_too_few_sequences.fasta', 'seq_type': 'DNA', 'error_text': LESS_THAN_TWO_SEQS_ERROR},
        {'file': 'protein_too_few_sequences.fasta', 'seq_type': 'Protein', 'error_text': LESS_THAN_TWO_SEQS_ERROR},
    ]

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
        self.assertTemplateUsed(response, 'base/index.html',
                                msg_prefix='Seq type: {}, file used: {}'.format(seq_type, input_file))

    def invalid_input_errors_are_shown_on_home_page(self, error_text, input_file='', seq_type='Protein'):
        """
        Tests that the correct error message is displayed on the home page on submission of an invalid alignment
        :param error_text: expected error message
        :param input_file: file containing invalid alignment
        :return:
        """
        response = self.response_for_invalid_post_request(input_file, seq_type)
        self.assertContains(response, escape(error_text),
                            msg_prefix='Seq type: {}, file used: {}'.format(seq_type, input_file))

    def invalid_input_passes_form_to_template(self, input_file='', seq_type='Protein'):
        """
        Tests that the form context is passed in response to an invalid alignment submission
        :param input_file: file containing invalid alignment
        :return:
        """
        response = self.response_for_invalid_post_request(input_file, seq_type)
        self.assertIsInstance(response.context['form'], QueryForm,
                              'Seq type: {}, file used: {}'.format(seq_type, input_file))

    def test_invalid_input(self):
        """
        Tests that submitting an invalid alignment input renders the home page, that the correct error message is
        displayed on the home page and that the form context is passed in response to an invalid alignment submission
        :return:
        """
        for file in self.files:
            # Does not work with Jinja2
            # self.invalid_input_renders_index_template(file['file'], file['seq_type'])

            self.invalid_input_errors_are_shown_on_home_page(file['error_text'], file['file'], file['seq_type'])

            # Does not work with Jinja2
            # self.invalid_input_passes_form_to_template(file['file'], file['seq_type'])


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
        align = Alignment.objects.create_alignment(name, data)
        self.response = self.client.get('/align-display/' + str(align.id) + '/')

    def test_align_display_page_uses_align_display_seq_template(self):
        """
        Tests that align_display view returns a 200 response on a GET request and uses the correct template
        """
        self.assertEqual(self.response.status_code, 200)

        # Does not work with Jinja2
        # self.assertTemplateUsed(self.response, 'base/align_display.html')

        with self.assertHTML(self.response, 'title[id="head-title"]') as elems:
            self.assertEqual('Formalign.eu Alignment Display', elems[0].text)

    def test_align_display_page_displays_default_formatted_protein_alignment(self):
        """
        Tests that align_display displays an alignment in the default format on a valid GET request (line of 80
        characters in blocks of 10 characters)
        """
        with self.assertHTML(self.response, 'tr[class="al_ln"]') as elems:
            self.assertEqual(len(
                [e for e in elems[0].findall('td') if e.attrib['class'] not in ['block_sep', 'display_artifact']]
            ),
                81, elems[0].getchildren()[0].text
            )
            self.assertEqual(len([e for e in elems[0].findall('td') if e.attrib['class'] == 'block_sep']),
                             8, elems[0].getchildren()[0].text
                             )

    def test_align_display_page_displays_correct_protein_alignment_sequence(self):
        """
        Tests that align_display displays an alignment with correct sequences
        """
        expected_seqs = file_to_string('spa_protein_alignment.fasta')
        align_expected = io.StringIO(expected_seqs)
        alignment = parse_fasta_alignment(align_expected)

        # get displayed sequences
        with self.assertHTML(self.response, 'tr') as elems:
            seq_disp = []
            for els in elems:
                seq_disp_line = []
                for e in els.findall('td')[:-1]:
                    if e.attrib['class'] in ['residue S0', 'residue S1']:
                        seq_disp_line.append(e.text)
                if seq_disp_line:
                    seq_disp.append(seq_disp_line)

        # recompose sequences
        re_seqs = []
        cat_re_seq = []
        for i in range(0, len(alignment) + 1):
            for j in range(i, len(seq_disp), len(alignment) + 1):
                re_seq = [seq_disp[j] for j in range(i, len(seq_disp), len(alignment) + 1)]
                cat_re_seq = []
                for r in re_seq:
                    cat_re_seq.extend(r)
            re_seqs.append(cat_re_seq)

        # check sequences against original alignment
        for i, al in enumerate(alignment):
            al_li = list(al.seq)
            self.assertEqual(al_li, re_seqs[i], re_seqs[i])

    def test_align_display_page_displays_sequences_in_the_correct_order(self):
        """
        Tests that align_display displays the sequences in the correct order
        """
        expected_seqs = file_to_string('spa_protein_alignment.fasta')
        align_expected = io.StringIO(expected_seqs)
        alignment = parse_fasta_alignment(align_expected)
        ids = [al.id for al in alignment]
        ids.append('consensus 70%')
        with self.assertHTML(self.response, 'tr') as elems:
            ids_disp = []
            for els in elems:
                for e in els.findall('td')[:-1]:
                    if e.attrib['class'] == 'seq_id':
                        ids_disp.append(e.text)

        # blocks of sequence IDs
        bl = [ids_disp[i:i + len(ids)] for i in range(0, len(ids_disp), len(ids))]
        for b in bl:
            for i in range(len(ids)):
                self.assertEqual(ids[i], b[i], b[i])

    def test_align_display_page_displays_correct_consensus(self):
        """
        Tests that align_display displays the correct consensus sequence
        """
        expected_seqs = file_to_string('spa_protein_alignment.fasta')
        align_expected = io.StringIO(expected_seqs)
        alignment = parse_fasta_alignment(align_expected)
        alignment = consensus_add(alignment)

        # get displayed sequences
        with self.assertHTML(self.response, 'tr') as elems:
            seq_disp = []
            for els in elems:
                seq_disp_line = []
                for e in els.findall('td')[:-1]:
                    if e.attrib['class'] in ['residue S0', 'residue S1']:
                        seq_disp_line.append(e.text)
                if seq_disp_line:
                    seq_disp.append(seq_disp_line)

        # recompose sequences
        cat_re_seq = []
        for j in range(len(alignment) - 1, len(seq_disp), len(alignment)):
            re_seq = [seq_disp[j] for j in range(len(alignment) - 1, len(seq_disp), len(alignment))]
            cat_re_seq = []
            for r in re_seq:
                cat_re_seq.extend(r)

        # check consensus
        cons_li = list(alignment[-1].seq)
        self.assertEqual(cons_li, cat_re_seq, cat_re_seq)

    def test_align_display_renders_correct_color_classes(self):
        """
        Tests that align_display assigns the correct color classes (residue S0 or residue S1) to the residues
        """
        expected_seqs = file_to_string('spa_protein_alignment.fasta')
        align_expected = io.StringIO(expected_seqs)
        alignment = parse_fasta_alignment(align_expected)
        alignment = consensus_add(alignment)
        alignment = annotate(alignment)

        # get displayed sequences
        with self.assertHTML(self.response, 'tr') as elems:
            seq_disp = []
            for els in elems:
                seq_disp_line = []
                for e in els.findall('td')[:-1]:
                    if e.attrib['class'] in ['residue S0', 'residue S1']:
                        seq_disp_line.append(e.attrib['class'])
                if seq_disp_line:
                    seq_disp.append(seq_disp_line)

        # recompose sequences
        re_seqs = []
        cat_re_seq = []
        for i in range(0, len(alignment)):
            for j in range(i, len(seq_disp), len(alignment)):
                re_seq = [seq_disp[j] for j in range(i, len(seq_disp), len(alignment))]
                cat_re_seq = []
                for r in re_seq:
                    cat_re_seq.extend(r)
            re_seqs.append(cat_re_seq)

        # check color classes
        for i, al in enumerate(alignment):
            al_li = ['residue S%s' % a for a in al.letter_annotations['eq']]
            self.assertEqual(al_li, re_seqs[i], re_seqs)


class AlignDisplayTestCaseSpeed(TestCase, AssertHTMLMixin):
    """
    Tests for alignment display
    """

    def setUp(self):
        """
        Creates a response from a GET request to /align-display/ with an alignment pk
        :param input_file: file containing alignment
        :return: response
        """
        name = 'SPA1 protein alignment'
        align_input = io.StringIO(file_to_string('spa1_protein_alignment.fasta'))
        data = parse_fasta_alignment(align_input)
        for d in data:
            d.seq.alphabet = Gapped(ExtendedIUPACProtein())
        self.align = Alignment.objects.create_alignment(name, data)

    def test_align_display_page_displays_protein_alignment_sequence_at_reasonable_speed(self):
        """
        Tests that align_display displays an alignment with correct sequences
        :return:
        """
        t0 = time.time()
        self.client.get('/align-display/' + str(self.align.id) + '/')
        resp_time = time.time() - t0
        self.assertTrue(resp_time < 2, 'response time: ' + format(resp_time))


class FormatHelpersTestCase(TestCase):
    """
    Tests for seq-display and align_display views helper functions
    """

    def setUp(self):
        """
        Creates an alignment from ser_thr_kin_short in the db
        """
        align_input = io.StringIO(file_to_string('ser_thr_kin_short.fasta'))
        data = parse_fasta_alignment(align_input)
        for d in data:
            d.seq.alphabet = Gapped(ExtendedIUPACProtein())
        self.alignment = data
        align_input_a = io.StringIO(file_to_string('protein_annotate_test.fasta'))
        data_a = parse_fasta_alignment(align_input_a)
        for d in data_a:
            d.seq.alphabet = Gapped(ExtendedIUPACProtein())
        self.alignment_a = data_a

    def test_consensus_get_returns_consensus_to_alignment(self):
        """
        Tests that consensus_add function returns correct default consensus
        """
        alignment = self.alignment
        cons_seq = AlignInfo.SummaryInfo(alignment).gap_consensus()
        cons_got = consensus_add(alignment)[-1]
        self.assertEqual(cons_seq, cons_got.seq, cons_got.seq)
        self.assertEqual('consensus 70%', cons_got.id, cons_got.id)

    def test_split_lines_alignment_splits_blocks_of_80_chars(self):
        """
        Tests that split_lines function splits blocks of alignment sequences in lines of 80 characters by default
        """
        alignment = self.alignment
        split_alignment = split_lines(alignment, split_type='alignment')
        for seq in split_alignment[:-1]:
            self.assertEqual([80 for i in range(len(seq))], [len(s) for s in seq])
        for seq in split_alignment[-1]:
            self.assertTrue(len(seq) <= 80)

    def test_split_lines_sequence_splits_lines_of_80_chars(self):
        """
        Tests that split_lines function splits sequences in lines of 80 characters by default
        """
        alignment = self.alignment
        split_sequence = split_lines(alignment, split_type='sequence')
        for seq in split_sequence:
            self.assertEqual([80 for i in range(len(seq['seq']) - 1)], [len(s) for s in seq['seq'][:-1]])
            self.assertTrue(len(seq['seq'][-1]) <= 80)

    def test_split_lines_in_blocks_returns_correct_object(self):
        """
        Tests that lines are split in blocks of 3 characters, tests that split_lines_in_blocks returns the correct
        object
        """
        split_seq_file = os.path.join(
            os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'test_data', 'split_seq.pckl'
        )
        split_sequence = pickle.load(open(split_seq_file, 'rb'))
        split_blocks = split_lines_in_blocks(split_sequence, block_size=3)
        split_blocks_li_exp = [
            [
                ['seq1', [[(0, 'C'), (1, 'N'), (0, 'G')], [(1, 'R'), (0, 'S'), (1, 'S')]]],
                ['seq2', [[(0, 'E'), (1, 'N'), (0, 'G')], [(1, 'R'), (0, 'G'), (1, 'S')]]],
                ['seq3', [[(0, 'D'), (0, 'R'), (0, 'D')], [(1, 'R'), (0, 'C'), (1, 'S')]]],
                ['seq4', [[(0, 'R'), (1, 'N'), (0, '-')], [(0, '-'), (0, '-'), (1, 'S')]]],
                ['consensus 70%', [[(0, 'X'), (1, 'N'), (0, 'X')], [(1, 'R'), (0, 'X'), (1, 'S')]]]
            ],
            [
                ['seq1', [[(0, '-'), (0, '-'), (0, '-')], [(0, '-')]]],
                ['seq2', [[(0, 'M'), (0, 'I'), (0, 'V')], [(0, 'S')]]],
                ['seq3', [[(0, '-'), (0, '-'), (0, '-')], [(0, '-')]]],
                ['seq4', [[(0, '-'), (0, '-'), (0, '-')], [(0, '-')]]],
                ['consensus 70%', [[(0, '-'), (0, '-'), (0, '-')], [(0, '-')]]]
            ]
        ]
        split_blocks_li = [
            [
                [line[0], [[(annot, res) for (annot, res) in res_block] for res_block in line[1]]
                 ] for line in seq_block
                ] for seq_block in split_blocks
            ]
        for seq_block in split_blocks_li:
            for line in seq_block:
                res_blocks = line[1]
                for res_block in res_blocks:
                    length = len([(annot, res) for (annot, res) in res_block])
                    self.assertTrue(length <= 3, length)
        self.assertEqual(split_blocks_li_exp, split_blocks_li, split_blocks_li)

    def test_annotate_returns_correct_default_consensus_annotation(self):
        """
        Tests that consensus is correctly annotated for display
        """
        exp_annot = [0, 1, 0, 1, 0, 1, 0, 0, 0, 0]
        alignment = consensus_add(self.alignment_a)
        annot_alignment = annotate(alignment)
        self.assertEqual(
            exp_annot,
            annot_alignment[-1].letter_annotations['eq'],
            annot_alignment[-1].letter_annotations['eq']
        )

    def test_annotate_returns_correct_sequence_annotations(self):
        """
        Tests that sequences are correctly annotated for display
        """
        exp_annot = [
            [0, 1, 0, 1, 0, 1, 0, 0, 0, 0],
            [0, 1, 0, 1, 0, 1, 0, 0, 0, 0],
            [0, 0, 0, 1, 0, 1, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 1, 0, 0, 0, 0],
            [0, 1, 0, 1, 0, 1, 0, 0, 0, 0]
        ]
        alignment = consensus_add(self.alignment_a)
        annot_alignment = annotate(alignment)
        for i, e in enumerate(exp_annot):
            self.assertEqual(
                e,
                annot_alignment[i].letter_annotations['eq'],
                annot_alignment[i].letter_annotations['eq']
            )
