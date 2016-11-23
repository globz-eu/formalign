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

from behave import when, then, use_step_matcher
import pyperclip
import re
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as ec
from helper_funcs.helpers_test import file_to_string
from lxml import html
from io import StringIO
from formalign.settings import TEST

__author__ = 'Stefan Dieterle'


use_step_matcher('re')


@when(
    r'the user pastes a (?P<sequence_type>.+) '
    r'alignment: "(?P<alignment_name>[^"]*)" '
    r'in the form text area'
)
def paste_alignment(context, sequence_type, alignment_name):
    """
    pastes an alignment in the form text area
    :param context: behave context
    :param sequence_type: sequence type
    :param alignment_name: alignment to submit
    """
    alignments = {
        'spa protein alignment': 'spa_protein_alignment.fasta',
        'empty': '',
        'invalid characters': sequence_type + '_invalid_characters.fasta',
        'too few sequences': sequence_type + '_too_few_sequences.fasta',
        'different sequence lengths': sequence_type + '_invalid_alignment.fasta',
        'invalid FASTA format': sequence_type + '_invalid_fasta.fasta',
    }
    alignment_string = file_to_string(alignments[alignment_name]) if alignments[alignment_name] else ''
    # pyperclip.copy(alignment_string)
    alignment_input = context.browser.find_element_by_css_selector('textarea#id_align_input')
    alignment_input.send_keys(alignment_string)
    # alignment_input.send_keys(Keys.CONTROL, 'v')


@when(r'a (?P<sequence_type>.*) alignment: "(?P<alignment_name>[^"]*)" is submitted')
def submit_alignment(context, sequence_type, alignment_name):
    """
    performs the request associated with submitting the requested alignment
    :param context: behave context
    :param sequence_type: sequence type
    :param alignment_name: alignment to submit
    """
    csrftoken = context.client.cookies['csrftoken']
    alignments = {
        'spa protein alignment': 'spa_protein_alignment.fasta',
        'empty': None,
        'invalid characters': sequence_type + '_invalid_characters.fasta',
        'too few sequences': sequence_type + '_too_few_sequences.fasta',
        'different sequence lengths': sequence_type + '_invalid_alignment.fasta',
        'invalid FASTA format': sequence_type + '_invalid_fasta.fasta',
    }
    alignment_string = file_to_string(alignments[alignment_name]) if alignments[alignment_name] else None
    context.client.headers.update({'referer': context.r.url})
    context.r = context.client.post(context.r.url,
                                    data={'csrfmiddlewaretoken': csrftoken, 'seq_type': sequence_type,
                                          'cons_type': 'identity',
                                          'align_input': alignment_string, 'custom_data': 'custom'})
    context.display = html.parse(StringIO(context.r.text)).getroot()


@then(r'the user (?:is redirected to|stays on|is on) the "(?P<page>[^"]*)" page')
def check_redirection(context, page):
    """
    tests that the user is redirected to / stays on the correct page
    :param context: behave context
    :param page: expected redirection page
    """
    pages = {
        'sequence display': 'Formalign.eu Sequence Display',
        'alignment display': 'Formalign.eu Alignment Display',
        'home': 'Formalign.eu Home',
        '404': 'Formalign.eu Error 404',
    }
    title = ''
    if TEST == 'acceptance':
        try:
            WebDriverWait(context.browser, 10).until(
                ec.title_is(pages[page])
            )
        finally:
            title = context.browser.title
    elif TEST == 'functional':
        title = context.display.cssselect('title[id="head-title"]')[0].text_content()
    assert pages[page] == title, 'Expected: %s\nGot: %s' % (pages[page], title)


@then(r'there (?:are|is a) (?P<sequence_type>.*) sequences? displayed')
def check_sequences_presence(context, sequence_type):
    """
    tests that there are sequences displayed on the query sequences page
    :param sequence_type: sequence type to check
    :param context: behave context
    """
    if sequence_type in ['protein', 'demo']:
        if TEST == 'acceptance':
            body_elem = context.browser.find_element_by_css_selector('body')
            body = context.browser.execute_script('return arguments[0].innerHTML', body_elem)
            body_html = html.parse(StringIO(body)).getroot()
            context.sequence_lines = body_html.cssselect('p[class="query_seq_display"]')
        elif TEST == 'functional':
            context.sequence_lines = context.display.cssselect('p[class="query_seq_display"]')
        assert context.sequence_lines, 'Got %s' % context.sequence_lines
    elif sequence_type == 'consensus':
        if TEST == 'acceptance':
            context.consensus_seq = context.browser.find_elements_by_css_selector(
                '.query_seq'
            )[-1].find_elements_by_css_selector(
                '.query_seq_display'
            )
        elif TEST == 'functional':
            context.consensus_seq = context.display.cssselect(
                'div[class="query_seq bg-color-body"]'
            )[-1].cssselect(
                'p[class="query_seq_display"]'
            )
        assert context.consensus_seq is not None, 'Got %s' % context.consensus_seq


@then(r'the sequences are displayed in lines of 80 characters')
def check_sequence_line_lengths(context):
    """
    tests that sequences are displayed in lines of 80 characters
    :param context: behave context
    """
    for f in context.sequence_lines:
        line_text = f.text
        assert len(line_text) <= 80, '%s was longer than 80 characters' % line_text


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
            line_text = context.sequence_lines[i].text_content()
            assert a == line_text,\
                'Sequence %s is not as expected' % line_text
        elif alignment_type == 'consensus':
            if TEST == 'acceptance':
                cons_line_text = context.consensus_seq[i].text
            elif TEST == 'functional':
                cons_line_text = context.consensus_seq[i].text_content()
            assert a == cons_line_text,\
                'Sequence %s is not as expected' % cons_line_text


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
    seqs_meta = []
    if TEST == 'acceptance':
        seqs_meta = context.browser.find_elements_by_css_selector('.query_seq_meta')
    elif TEST == 'functional':
        seqs_meta = context.display.find_class('query_seq_meta')
    meta_text = ''
    if alignment_type in ['demo', 'protein']:
        for i, a in enumerate(seqs_meta_expected):
            if TEST == 'acceptance':
                meta_text = seqs_meta[i].text
            elif TEST == 'functional':
                meta_text = seqs_meta[i].text_content()
            assert a == meta_text, 'Metadata %s is not expected' % meta_text
    elif alignment_type == 'consensus':
        if TEST == 'acceptance':
            meta_text = seqs_meta[-1].text
        elif TEST == 'functional':
            meta_text = seqs_meta[-1].text_content()
        assert seqs_meta_expected[0] == meta_text,\
            'Expected: %s\nGot: %s' % (seqs_meta_expected[0], meta_text)


@then(r'there is a "(?P<button>[^"]*)" button with "(?P<method>[^"]*)" method and "(?P<action>[^"]*)" action')
def check_render_button(context, button, method, action):
    """
    tests that the button has the expected attributes
    :param context: behave context
    :param button: button name
    :param method: action method
    :param action: action URL
    """
    render_form_method = ''
    if TEST == 'acceptance':
        render_form_method = context.browser.find_element_by_id('render').get_attribute('method')
        context.render_form_action = context.browser.find_element_by_id('render').get_attribute('action')
    if TEST == 'functional':
        render_form_method = context.display.cssselect('form[id="render"]')[0].attrib.get('method')
        context.render_form_action = context.display.cssselect('form[id="render"]')[0].attrib.get('action')
    assert method == render_form_method, 'Expected: %s\nGot: %s' % (method, render_form_method)
    assert action == context.render_form_action.split('/')[-3], \
        'Expected: %s\nGot: %s' % (action, context.render_form_action.split('/')[-3])


@then(r'the action URL of the "(?P<button>[^"]*)" button contains a 16 character slug')
def check_slug_in_render_button_action(context, button):
    """
    tests that action URL in button contains a 16 character slug
    :param context: behave context
    :param button: button to check
    """
    slug_pattern = re.compile('^([a-zA-Z]|\d){16}$')
    assert re.match(slug_pattern, context.render_form_action.split('/')[-2]), \
        '%s did not match slug pattern' % context.render_form_action.split('/')[-2]
