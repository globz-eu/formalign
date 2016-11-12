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

from behave import when, then
import pyperclip
from selenium.webdriver.common.keys import Keys
from helper_funcs.helpers_test import file_to_string

__author__ = 'Stefan Dieterle'


@when(r'the user pastes a (?P<alignment_type>.+) alignment in the form text area')
def paste_alignment(context, alignment_type):
    """
    pastes an alignment in the form text area
    :param context: behave context
    :param alignment_type: alignment to paste in
    """
    alignments = {'protein': 'spa_protein_alignment.fasta'}
    alignment_string = file_to_string(alignments[alignment_type])
    pyperclip.copy(alignment_string)
    alignment_input = context.browser.find_element_by_css_selector('textarea#id_align_input')
    alignment_input.send_keys(Keys.CONTROL, 'v')


@then(r'the user is redirected to the "(?P<page>[^"]*)" page')
def check_redirection(context, page):
    """
    tests that the user is redirected to the correct page
    :param context: behave context
    :param page: expected redirection page
    """
    pages = {
        'sequence display': 'Formalign.eu Sequence Display',
        'alignment display': 'Formalign.eu Alignment Display'
    }
    assert pages[page] == context.browser.title, 'Got %s' % context.browser.title


@then(r'there (?:are|is a) (?P<sequence_type>.*) sequences? displayed')
def check_sequences_presence(context, sequence_type):
    """
    tests that there are sequences displayed
    :param sequence_type: sequence type to check
    :param context: behave context
    """
    if sequence_type in ['protein', 'demo']:
        context.sequence_lines = context.browser.find_elements_by_css_selector('.query_seq_display')
        assert context.sequence_lines, 'Got %s' % context.sequence_lines
    elif sequence_type == 'consensus':
        context.consensus_seq = context.browser.find_elements_by_css_selector('.query_seq')[-1].find_elements_by_css_selector('.query_seq_display')
        assert context.consensus_seq is not None, 'Got %s' % context.consensus_seq


@then(r'the sequences are displayed in lines of 80 characters')
def check_sequence_line_lengths(context):
    """
    tests that sequences are displayed in lines of 80 characters
    :param context: behave context
    """
    for f in context.sequence_lines:
        assert len(f.text) <= 80, '%s was longer than 80 characters' % f.text


@then(r'the correct (?P<alignment_type>.+) sequences? (?:are|is) displayed')
def check_correct_sequences(context, alignment_type):
    """
    tests that the sequence display page displays the correct sequences
    :param context: behave context
    :param alignment_type: alignment type for determining correct first sequence
    """
    displayed_seqs = {
        'demo': 'ser_thr_kinase_family_display.txt',
        'protein': 'spa_protein_alignment_display.txt',
        'consensus': 'consensus_display.txt',
    }
    seqs = file_to_string(displayed_seqs[alignment_type]).splitlines()
    for i, a in enumerate(seqs):
        if alignment_type in ['demo', 'protein']:
            assert a == context.sequence_lines[i].text, 'Sequence %s is not as expected' % context.sequence_lines[i].text
        elif alignment_type == 'consensus':
            assert a == context.consensus_seq[i].text, 'Sequence %s is not as expected' % context.consensus_seq[i].text


@then(r'the correct (?P<alignment_type>.+) sequence metadata (?:are|is) displayed')
def check_correct_sequence_metadata(context, alignment_type):
    """
    tests that the sequence display page displays the correct sequence metadata
    :param context: behave context
    :param alignment_type: alignment type for determining correct first sequence
    """
    displayed_seqs_meta = {
        'demo': 'ser_thr_kinase_family_display_meta.txt',
        'protein': 'spa_protein_alignment_display_meta.txt',
        'consensus': 'consensus_display_meta.txt',
    }
    seqs_meta_expected = file_to_string(displayed_seqs_meta[alignment_type]).splitlines()
    seqs_meta = context.browser.find_elements_by_css_selector('.query_seq_meta')
    if alignment_type in ['demo', 'protein']:
        for i, a in enumerate(seqs_meta_expected):
            assert a == seqs_meta[i].text, 'Metadata %s is not expected' % seqs_meta[i].text
    elif alignment_type == 'consensus':
        assert seqs_meta_expected[0] == seqs_meta[-1].text,\
            'Got %s, expected %s' % (seqs_meta[-1].text, seqs_meta_expected[0])
