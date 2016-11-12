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

import re

from behave import when, then

from lxml import html
from io import StringIO

from helper_funcs.helpers_test import file_to_string

__author__ = 'Stefan Dieterle'


@when(r'the "(?P<button>[^"]*)" button is pressed')
def click_button(context, button):
    """
    performs the request associated with clicking the reuqested button
    :param context: behave context
    :param button: requested button
    """
    if button == 'Demo':
        csrftoken = context.client.cookies['csrftoken']
        context.client.headers.update({'referer': context.home_url + '/'})
        context.r = context.client.post(context.home_url + '/',
                                        data={'csrfmiddlewaretoken': csrftoken, 'seq_type': 'DNA',
                                              'align_input': '', 'custom_data': 'demo'})
    elif button == 'Render':
        context.r = context.client.get(
            context.home_url + context.display.cssselect('form[id="render"]')[0].attrib.get('action')
        )
    context.display = html.parse(StringIO(context.r.text)).getroot()


@when(r'a custom (?P<sequence_type>.*) alignment: "(?P<alignment_name>[^"]*)" is submitted')
def submit_custom_alignment(context, sequence_type, alignment_name):
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


@then(r'the user is redirected to the "(?P<page>[^"]*)" page')
def check_redirect(context, page):
    """
    tests that the user is redirected to the correct page
    :param context: behave context
    :param page: expected redirection page
    """
    title = context.display.cssselect('title[id="head-title"]')
    pages = {
        'sequence display': {
            'title': 'Formalign.eu Sequence Display',
            'url': '%s/query-sequences'
        },
        'alignment display': {
            'title': 'Formalign.eu Alignment Display',
            'url': '%s/align-display'
        }
    }
    assert pages[page]['title'] == title[0].text_content(), 'Got %s' % title[0].text_content()
    assert pages[page]['url'] % context.home_url == '/'.join(context.r.url.split('/')[:-2])


@then(r'there (?:are|is a) (?P<sequence_type>.*) sequences? displayed')
def check_sequences_presence(context, sequence_type):
    """
    tests that there are sequences displayed
    :param sequence_type: sequence type to check
    :param context: behave context
    """
    if sequence_type in ['protein', 'demo']:
        context.sequence_lines = context.display.cssselect('p[class="query_seq_display"]')
        assert context.sequence_lines, 'Got %s' % context.sequence_lines
    elif sequence_type == 'consensus':
        context.consensus_seq = context.display.cssselect(
            'div[class="query_seq bg-color-body"]'
        )[-1].cssselect(
            'p[class="query_seq_display"]'
        )
        assert context.consensus_seq is not None, 'Got %s' % context.consensus_seq


@then(r'the sequences are displayed in lines of 80 characters')
def check_sequences_line_length(context):
    """
    tests that sequences are displayed in lines of 80 characters
    :param context: behave context
    """
    for l in context.sequence_lines:
        assert len(l.text_content()) <= 80, 'Line %s has length %d' % (l.text_content(), len(l.text_content()))


@then(r'the correct (?P<alignment_type>.*) sequences? (?:are|is) displayed')
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
            assert a == context.sequence_lines[i].text_content(), \
                'Sequence %s is not as expected' % context.sequence_lines[i].text_content()
        elif alignment_type == 'consensus':
            assert a == context.consensus_seq[i].text_content(), \
                'Sequence %s is not as expected' % context.consensus_seq[i].text_content()


@then(r'the correct (?P<alignment_type>.*) sequence metadata (?:are|is) displayed')
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
    seqs_meta = context.display.find_class('query_seq_meta')
    if alignment_type in ['demo', 'protein']:
        for i, a in enumerate(seqs_meta_expected):
            assert a == seqs_meta[i].text_content(), 'Metadata %s is not expected' % seqs_meta[i].text_content()
    elif alignment_type == 'consensus':
        assert seqs_meta_expected[0] == seqs_meta[-1].text_content(),\
            'Got %s, expected %s' % (seqs_meta[-1].text_content(), seqs_meta_expected[0])


@then(r'there is a "(?P<button>[^"]*)" button with "(?P<method>[^"]*)" method and "(?P<action>[^"]*)" action')
def check_render_button(context, button, method, action):
    """
    tests that the button has the expected attributes
    :param context: behave context
    :param button: button name
    :param method: action method
    :param action: action URL
    """
    context.render_form = context.display.cssselect('form[id="render"]')
    assert method == context.render_form[0].attrib.get('method'), 'Got %s' % context.render_form[0].attrib.get('method')
    assert action == context.render_form[0].attrib.get('action').split('/')[1], \
        'Got %s' % context.render_form[0].attrib.get('action').split('/')[1]


@then(r'the action URL of the "(?P<button>[^"]*)" button contains a 16 character slug')
def check_slug_in_render_button_action(context, button):
    """
    tests that action URL in button contains a 16 character slug
    :param context: behave context
    :param button: button to check
    """
    slug_pattern = re.compile('^([a-zA-Z]|\d){16}$')
    assert re.match(slug_pattern, context.render_form[0].attrib.get('action').split('/')[-2]), \
        '%s did not match slug pattern' % context.render_form[0].attrib.get('action').split('/')[-2]
