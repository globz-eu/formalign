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

from behave import then, use_step_matcher
from formalign.settings import TEST
from helper_funcs.helpers_test import file_to_string
from lxml import html
from io import StringIO
from functional_tests.test_implementation.alignment_rendering import alignment_formatting, get_displayed_seqs

__author__ = 'Stefan Dieterle'


use_step_matcher('re')


@then(r'the alignment is displayed with 80 characters per line in blocks of 10 with sequence IDs')
def check_alignment_formatting(context):
    """
    tests that the alignments are displayed with the correct formatting
    :param context: behave context
    """
    tables = None
    seqs_meta = file_to_string('spa_protein_alignment_meta.txt').splitlines()
    if TEST == 'acceptance':
        align_display_html = context.browser.find_element_by_css_selector('.align_display').get_attribute('innerHTML')
        align_display = html.parse(StringIO(align_display_html)).getroot()
        tables = align_display.find_class('align_table')
    elif TEST == 'functional':
        tables = context.display.find_class('align_table')

    alignment_formatting(seqs_meta, tables)


@then(r'the expected alignments are displayed')
def check_alignment_sequences(context):
    """
    tests that the expected sequences are displayed
    :param context: behave context
    """
    alignment = file_to_string('spa_protein_alignment_seqs.txt')
    alignment_list = [[c for c in a] for a in alignment.split('\n')[:-1]]
    elems = None
    if TEST == 'acceptance':
        align_display = context.browser.find_element_by_css_selector('.align_display').get_attribute('innerHTML')
        display = html.parse(StringIO(align_display)).getroot()
        elems = display.cssselect('tr')
    elif TEST == 'functional':
        elems = context.display.cssselect('tr')

    re_seqs = get_displayed_seqs(elems, len(alignment_list))

    for i, al_li in enumerate(alignment_list):
        assert al_li == re_seqs[i], 'expected: %s\n got: %s' % (al_li, re_seqs[i])


@then(r'the expected consensus sequence is displayed')
def check_consensus_sequence(context):
    """
    tests that the expected consensus sequence is displayed
    :param context: behave context
    """
    alignment = file_to_string('spa_protein_alignment_consens.txt')
    alignment_list = [[c for c in a] for a in alignment.split('\n')[:-1]]
    elems = None
    if TEST == 'acceptance':
        align_display = context.browser.find_element_by_css_selector('.align_display').get_attribute('innerHTML')
        display = html.parse(StringIO(align_display)).getroot()
        elems = display.cssselect('tr')
    elif TEST == 'functional':
        elems = context.display.cssselect('tr')

    cat_re_seq = get_displayed_seqs(elems, len(alignment_list), cons=True)

    cons_li = alignment_list[-1]
    assert cons_li == cat_re_seq, cat_re_seq


@then(r'the sequence elements have the expected color classes')
def check_alignment_sequences_annotation(context):
    """
    tests that the sequence elements (residues or bases) have the expected
    color classes
    :param context: behave context
    """
    alignment = file_to_string('spa_protein_alignment_seqs_annot.txt')
    alignment_list = [['residue S%s' % a for a in al] for al in alignment.split('\n')[:-1]]
    elems = None
    if TEST == 'acceptance':
        align_display = context.browser.find_element_by_css_selector('.align_display').get_attribute('innerHTML')
        display = html.parse(StringIO(align_display)).getroot()
        elems = display.cssselect('tr')
    elif TEST == 'functional':
        elems = context.display.cssselect('tr')

    re_seqs = get_displayed_seqs(elems, len(alignment_list), annot=True)

    for i, al in enumerate(alignment_list):
        assert al == re_seqs[i], 'expected: %s\n got: %s' % (al, re_seqs[i])
