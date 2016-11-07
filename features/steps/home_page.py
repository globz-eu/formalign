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

from behave import given, when, then

from formalign.settings import SERVER_URL
from lxml import html
from io import StringIO

__author__ = 'Stefan Dieterle'


@given(r'a user visits the URL "(?P<url>[^"]*)"')
def visit_url(context, url):
    if SERVER_URL == 'liveserver':
        context.home_url = context.base_url
    else:
        context.home_url = SERVER_URL
    context.r = context.client.get(context.home_url + url)


@when(r'the user looks at the page')
def get_elements(context):
    context.display = html.parse(StringIO(context.r.text)).getroot()


@then(r'the server\'s response status code is (?P<expected_code>\d+)')
def check_status_code(context, expected_code):
    assert context.r.status_code == int(expected_code), 'Got %d' % context.r.status_code


@then(r'the url is the home url')
def check_page_url(context):
    assert context.r.url == context.home_url + '/', 'Got %s' % context.r.url


@then(r'the page title is "(?P<expected_title>[^"]*)"')
def check_page_title(context, expected_title):
    title = context.display.cssselect('title[id="head-title"]')[0].text_content()
    assert expected_title == title, 'Got %s' % title


@then(r'the brand text says "(?P<expected_brand>[^"]*)"')
def check_brand(context, expected_brand):
    brand = context.display.cssselect('a[class="navbar-brand"]')[0].text_content()
    assert expected_brand == brand, 'Got %s' % brand


@then(r'there is a form labeled "(?P<expected_form_label>[^"]*)"')
def check_form_label(context, expected_form_label):
    form_label = context.display.cssselect('textarea[id="id_align_input"]')[0].label.text_content()
    assert expected_form_label == form_label, 'Got %s' % form_label


@then(r'there is a text area with a placeholder saying "(?P<expected_textarea_placeholder>[^"]*)"')
def check_form_textarea_placeholder(context, expected_textarea_placeholder):
    textarea_placeholder = context.display.cssselect('textarea[id="id_align_input"]')[0].attrib.get('placeholder')
    assert expected_textarea_placeholder == textarea_placeholder, 'Got %s' % textarea_placeholder


@then(r'there are radio buttons labeled "(?P<expected_sequence_type_label>[^"]*)"')
def check_sequence_type_label(context, expected_sequence_type_label):
    sequence_type_label = context.display.cssselect('input[id="id_seq_type_0"]')[0].label.text_content()
    assert expected_sequence_type_label == sequence_type_label, 'Got %s' % sequence_type_label


@then(r'there is a "(?P<expected_radio_button_text>[^"]*)" (?P<choice>[^"]*) type radio button')
def check_sequence_type_radio_buttons(context, expected_radio_button_text, choice):
    radio_button = {
        'input sequence': {
            'DNA': context.display.cssselect('input[id="id_seq_type_1"]')[0].label.text_content(),
            'Protein': context.display.cssselect('label[for="id_seq_type_0"]')[1].text_content()
        },
        'consensus': {
            '% identity': context.display.cssselect('input[id="id_cons_type_1"]')[0].label.text_content(),
            'substitution matrix': context.display.cssselect('label[for="id_cons_type_0"]')[1].text_content()
        },
    }
    radio_button_text = radio_button[choice][expected_radio_button_text]
    assert re.match(
        '^\s+' + re.escape(expected_radio_button_text) + '$',
        radio_button_text
    ), 'Got %s' % radio_button_text


@then(r'the "(?P<type>[^"]*)" button is checkable')
def check_type_buttons_are_checkable(context, type):
    radio_button = {
        'DNA': '1',
        'Protein': '0'
    }
    button_checkable = context.display.cssselect('input[id="id_seq_type_%s"]' % radio_button[type])[0].checkable
    assert button_checkable, 'Got %s' % button_checkable


@then(r'the "(?P<type>[^"]*)" button is(?P<default>.*) checked by default')
def check_type_buttons_are_in_the_right_default_state(context, type, default):
    radio_button = {
        'DNA': '1',
        'Protein': '0'
    }
    if default == '':
        test = True
    elif default == ' not':
        test = False
    else:
        test = 'some unlikely string'
    button_checked = context.display.cssselect('input[id="id_seq_type_%s"]' % radio_button[type])[0].checked
    assert button_checked == test


@then(
    r'there is a "(?P<button_text>[^"]*)"'
    r' button named "(?P<button_name>[^"]*)"'
    r' with the value "(?P<button_value>[^"]*)"'
)
def check_button_text(context, button_text, button_name, button_value):
    button_id = {'Demo': 'submit-demo', 'Submit': 'submit-align'}
    button = context.display.cssselect('button[id="%s"]' % button_id[button_text])
    assert re.match(
        '^\\n\s+' + re.escape(button_text) + '\\n\s+$', button[0].text_content()
    ), 'Got %s' % button[0].text_content()
    assert button[0].attrib.get('name') == button_name, 'Got %s' % button[0].attrib.get('name')
    assert button[0].attrib.get('value') == button_value, 'Got %s' % button[0].attrib.get('value')


@then(r'there is a "(?P<link_text>[^"]*)" button with "(?P<href>[^"]*)" href')
def check_link_href(context, link_text, href):
    home_button = context.display.cssselect('a[class="navbar-brand"]')
    assert link_text == home_button[0].text_content(), 'Got %s' % home_button[0].text_content()
    assert href == home_button[0].attrib.get('href'), 'Got %s' % home_button[0].attrib.get('href')
