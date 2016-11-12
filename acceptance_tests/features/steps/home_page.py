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


__author__ = 'Stefan Dieterle'


@given(r'a user visits the URL "(?P<url>[^"]*)"')
def visit_url(context, url):
    """
    browses to given URL
    :param context: behave context
    :param url: URL to visit
    """
    if SERVER_URL == 'liveserver':
        context.home_url = context.base_url
    else:
        context.home_url = SERVER_URL
    context.browser.get(context.home_url + url)


@when(r'the user clicks the "(?P<choice_type>[^"]*)" (?P<button_type>.+) button')
def click_button(context, choice_type, button_type):
    """
    clicks the requested button
    :param context: behave context
    :param choice_type: button name
    :param button_type: button type
    """
    radio_button = {
        'DNA': {'type': 'seq', 'number': '1'},
        'Protein': {'type': 'seq', 'number': '0'},
        'Identity': {'type': 'cons', 'number': '0'},
        'Substitution Matrix': {'type': 'cons', 'number': '1'}
    }
    submit_button = {'Submit': 'submit-align', 'Render': 'render-align'}
    if button_type == 'radio':
        button = context.browser.find_element_by_css_selector(
            'input[id="id_%s_type_%s"]' % (
                radio_button[choice_type]['type'],
                radio_button[choice_type]['number']
            )
        )
    elif button_type == 'submit':
        button = context.browser.find_element_by_id(submit_button[choice_type])
    button.click()


# TODO: adapt to other URLs and implement on other pages
@then(r'the url is the home url')
def check_page_url(context):
    """
    tests that the current URL is the home URL
    :param context: behave context
    """
    assert context.browser.current_url == context.home_url + '/', 'Got %s' % context.r.url


@then(r'the page title is "(?P<expected_title>[^"]*)"')
def check_page_title(context, expected_title):
    """
    tests that the page title is the expected title
    :param context: behave context
    :param expected_title: expected title
    """
    assert expected_title == context.browser.title, 'Got %s' % context.browser.title


@then(r'the brand text says "(?P<expected_brand>[^"]*)"')
def check_brand(context, expected_brand):
    """
    tests that the brand text is as expected
    :param context: behave context
    :param expected_brand: expected brand text
    """
    brand = context.browser.find_element_by_css_selector('.navbar-brand').text
    assert expected_brand == brand, 'Got %s' % brand


@then(r'there is a form with a label matching "(?P<expected_form_label>[^"]*)"')
def check_form_label(context, expected_form_label):
    """
    tests that the form has the expected label
    :param context: behave context
    :param expected_form_label: expected label text
    """
    form_label = context.browser.find_element_by_css_selector('label[for="id_align_input"]').text
    assert re.match(expected_form_label, form_label), 'Got %s' % form_label


@then(r'there is a text area with a placeholder saying "(?P<expected_textarea_placeholder>[^"]*)"')
def check_form_textarea_placeholder(context, expected_textarea_placeholder):
    """
    tests that the form text area has the expected placeholder
    :param context: behave context
    :param expected_textarea_placeholder: expected text area placeholder text
    """
    textarea_placeholder = context.browser.find_element_by_css_selector(
        'textarea#id_align_input'
    ).get_attribute('placeholder')
    assert expected_textarea_placeholder == textarea_placeholder, 'Got %s' % textarea_placeholder


@then(r'there are radio buttons labeled "(?P<expected_sequence_type_label>[^"]*)"')
def check_sequence_type_label(context, expected_sequence_type_label):
    """
    tests the presence of radio button
    :param context: behave context
    :param expected_sequence_type_label: expected label
    """
    sequence_type_label = context.browser.find_elements_by_css_selector('label[for="id_seq_type_0"]')[0].text
    assert expected_sequence_type_label == sequence_type_label, 'Got %s' % sequence_type_label


@then(r'there is a "(?P<expected_radio_button_text>[^"]*)" (?P<choice>[^"]*) type radio button')
def check_sequence_type_radio_buttons(context, expected_radio_button_text, choice):
    """
    tests that the text of the radio button is as expected
    :param context: behave context
    :param expected_radio_button_text: expected radio button text
    :param choice: radio button type
    """
    radio_button = {
        'input sequence': {
            'DNA': context.browser.find_elements_by_css_selector('label[for="id_seq_type_1"]')[0].text,
            'Protein': context.browser.find_elements_by_css_selector('label[for="id_seq_type_0"]')[1].text
        },
        'consensus': {
            'Identity': context.browser.find_elements_by_css_selector('label[for="id_cons_type_0"]')[1].text,
            'Substitution Matrix': context.browser.find_elements_by_css_selector('label[for="id_cons_type_1"]')[0].text
        },
    }
    radio_button_text = radio_button[choice][expected_radio_button_text]
    assert expected_radio_button_text == radio_button_text, 'Got %s' % radio_button_text


@then(r'the "(?P<button_type>[^"]*)" button is(?P<default>.*) checked')
def check_type_buttons_are_in_the_right_default_state(context, button_type, default):
    """
    tests whether a radio button has the expected status (checked or not)
    :param context: behave context
    :param button_type: button to check
    :param default: button state
    """
    radio_button = {
        'DNA': {'type': 'seq', 'number': '1'},
        'Protein': {'type': 'seq', 'number': '0'},
        'Identity': {'type': 'cons', 'number': '0'},
        'Substitution Matrix': {'type': 'cons', 'number': '1'}
    }
    button_checked = context.browser.find_element_by_css_selector(
        'input[id="id_%s_type_%s"]' % (
            radio_button[button_type]['type'],
            radio_button[button_type]['number']
        )
    ).get_attribute('checked')
    if default == '':
        assert button_checked == 'true', 'Got %s' % button_checked
    elif default == ' not':
        assert not button_checked, 'Got %s' % button_checked


@then(
    r'there is a "(?P<button_text>[^"]*)"'
    r' button named "(?P<button_name>[^"]*)"'
    r' with the value "(?P<button_value>[^"]*)"'
)
def check_button_text(context, button_text, button_name, button_value):
    """
    tests presence and attributes of a button
    :param context: behave context
    :param button_text: button text
    :param button_name: button name
    :param button_value: button value
    """
    button_id = {'Demo': 'submit-demo', 'Submit': 'submit-align'}
    button = context.browser.find_element_by_css_selector('button[id="%s"]' % button_id[button_text])
    assert button_text == button.text, 'Got %s' % button.text
    assert button.get_attribute('name') == button_name, 'Got %s' % button.get_attribute('name')
    assert button.get_attribute('value') == button_value, 'Got %s' % button.get_attribute('value')


@then(r'there is a "(?P<link_text>[^"]*)" button with "(?P<href>[^"]*)" href')
def check_link_href(context, link_text, href):
    """
    tests the presence and attributes of a link
    :param context: behave context
    :param link_text: link text
    :param href: URL of link
    """
    home_button = context.browser.find_elements_by_css_selector('a[class="navbar-brand"]')
    assert link_text == home_button[0].text, 'Got %s' % home_button[0].text
    assert context.home_url + href == home_button[0].get_attribute('href'),\
        'Got %s' % home_button[0].get_attribute('href')
