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

from selenium import webdriver
from selenium.webdriver.common.desired_capabilities import DesiredCapabilities
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as ec
from formalign.settings import SERVER_URL, TEST, CHROME_DRIVER, FIREFOX_BINARY
import re
from behave import given, when, then, use_step_matcher
from lxml import html
from io import StringIO

__author__ = 'Stefan Dieterle'

use_step_matcher('re')


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
    context.r = context.client.get(context.home_url + url)


@given(r'a user visits the URL "(?P<url>[^"]*)" with "(?P<browser>[^"]*)"')
def visit_url_with_browser(context, url, browser):
    """
    browses to given URL
    :param context: behave context
    :param url: URL to visit
    :param browser: browser used to visit URL
    """
    context.browser_name = browser
    if browser == 'Chrome':
        if SERVER_URL == 'liveserver':
            context.browser = webdriver.Chrome(CHROME_DRIVER)
        else:
            caps = webdriver.DesiredCapabilities.CHROME
            context.browser = webdriver.Remote('http://selenium-hub:4444/wd/hub', caps.copy())
    elif browser == 'Firefox':
        if SERVER_URL == 'liveserver':
            caps = DesiredCapabilities.FIREFOX
            caps['marionette'] = True
            caps['binary'] = FIREFOX_BINARY
            context.browser = webdriver.Firefox(capabilities=caps)
        else:
            caps = webdriver.DesiredCapabilities.FIREFOX
            caps['marionette'] = True
            context.browser = webdriver.Remote('http://selenium-hub:4444/wd/hub', caps.copy())
    if SERVER_URL == 'liveserver':
        context.home_url = context.base_url
    else:
        context.home_url = SERVER_URL
    context.response = context.browser.get(context.home_url + url)


@given(r'active and inactive buttons')
def define_buttons_status(context):
    """
    get list of active and inactive buttons to be tested
    :param context: behave context
    """
    context.buttons_status = context.table


@when(r'the user looks at the page')
def get_elements(context):
    """
    parses the displayed page's HTML
    :param context: behave context
    """
    context.display = html.parse(StringIO(context.r.text)).getroot()


@then(r'the server\'s response status code is (?P<expected_code>\d+)')
def check_status_code(context, expected_code):
    """
    tests that the server response contains the expected status code
    :param context: behave context
    :param expected_code: expected status code
    """
    assert context.r.status_code == int(expected_code), 'Got %d' % context.r.status_code


@when(r'the user clicks the "(?P<choice_type>[^"]*)" (?P<button_type>.+) button')
def click_button(context, choice_type, button_type):
    """
    clicks the requested button or performs the request associated with
    clicking the requested button
    :param context: behave context
    :param choice_type: button name
    :param button_type: button type
    """
    radio_button = {
        'dna': {'type': 'seq', 'number': '1'},
        'protein': {'type': 'seq', 'number': '0'},
        'identity': {'type': 'cons', 'number': '0'},
        'substitution matrix': {'type': 'cons', 'number': '1'}
    }
    submit_button = {
        'submit': {'id': 'submit-align', 'redir_title': 'Formalign.eu Sequence Display'},
        'render': {'id': 'render-align', 'redir_title': 'Formalign.eu Alignment Display'},
        'demo': {'id': 'submit-demo', 'redir_title': 'Formalign.eu Sequence Display'}
    }
    if TEST == 'acceptance':
        if button_type == 'radio':
            button = context.browser.find_element_by_css_selector(
                'input[id="id_%s_type_%s"]' % (
                    radio_button[choice_type.lower()]['type'],
                    radio_button[choice_type.lower()]['number']
                )
            )
            button.click()
            try:
                WebDriverWait(context.browser, 10).until(
                    ec.element_to_be_selected(button)
                )
            finally:
                pass
        elif button_type == 'submit':
            button = context.browser.find_element_by_id(submit_button[choice_type.lower()]['id'])
            button.click()

    elif TEST == 'functional' and button_type == 'submit':
        if choice_type == 'Demo':
            csrftoken = context.client.cookies['csrftoken']
            context.client.headers.update({'referer': context.home_url + '/'})
            context.r = context.client.post(context.home_url + '/',
                                            data={'csrfmiddlewaretoken': csrftoken, 'seq_type': 'DNA',
                                                  'align_input': '', 'custom_data': 'demo'})
        if choice_type == 'Render':
            context.r = context.client.get(
                context.home_url + context.display.cssselect('form[id="render"]')[0].attrib.get('action')
            )
        context.display = html.parse(StringIO(context.r.text)).getroot()


@when(r'the user clicks the active radio buttons')
@when(r'the user looks at the radio buttons')
def define_radio_buttons(context):
    """
    defines active and inactive radio buttons in context
    :param context: behave context
    """
    radio_buttons = {
        'dna': {'type': 'seq', 'number': '1'},
        'protein': {'type': 'seq', 'number': '0'},
        'identity': {'type': 'cons', 'number': '0'},
        'substitution matrix': {'type': 'cons', 'number': '1'}
    }
    context.active_radio_buttons = [radio_buttons[row['active'].lower()] for row in context.buttons_status]
    context.inactive_radio_buttons = [radio_buttons[row['inactive'].lower()] for row in context.buttons_status]


@then(r'the current URL is the "(?P<expected_url>[^"]*)" URL')
def check_page_url(context, expected_url):
    """
    tests that the current URL is the home URL
    :param expected_url: expected current URL
    :param context: behave context
    """
    urls = {
        'home': context.home_url + '/',
        'sequence display': context.home_url + '/query-sequences/',
        'alignment display': context.home_url + '/align-display/'
    }
    if TEST == 'acceptance':
        assert urls[expected_url] == context.browser.current_url[:len(urls[expected_url])], \
            'Got %s' % context.browser.current_url
    elif TEST == 'functional':
        assert urls[expected_url] == context.r.url[:len(urls[expected_url])], \
            'Got %s' % context.r.url


@then(r'the brand text says "(?P<expected_brand>[^"]*)"')
def check_brand(context, expected_brand):
    """
    tests that the brand text is as expected
    :param context: behave context
    :param expected_brand: expected brand text
    """
    brand = ''
    if TEST == 'acceptance':
        brand = context.browser.find_element_by_css_selector('.navbar-brand').text
    elif TEST == 'functional':
        brand = context.display.cssselect('a[class="navbar-brand"]')[0].text_content()
    assert expected_brand == brand, 'Got %s' % brand


@then(r'there is a form with a label matching "(?P<expected_form_label>[^"]*)"')
def check_form_label(context, expected_form_label):
    """
    tests that the form has the expected label
    :param context: behave context
    :param expected_form_label: expected label text
    """
    form_label = ''
    if TEST == 'acceptance':
        form_label = context.browser.find_element_by_css_selector('label[for="id_align_input"]').text
    elif TEST == 'functional':
        form_label = context.display.cssselect('textarea[id="id_align_input"]')[0].label.text_content()
    assert re.match(expected_form_label, form_label), 'Got %s' % form_label


@then(r'there is a text area with a placeholder saying "(?P<expected_textarea_placeholder>[^"]*)"')
def check_form_textarea_placeholder(context, expected_textarea_placeholder):
    """
    tests that the form text area has the expected placeholder
    :param context: behave context
    :param expected_textarea_placeholder: expected text area placeholder text
    """
    textarea_placeholder = ''
    if TEST == 'acceptance':
        textarea_placeholder = context.browser.find_element_by_css_selector(
            'textarea#id_align_input'
        ).get_attribute('placeholder')
    elif TEST == 'functional':
        textarea_placeholder = context.display.cssselect('textarea[id="id_align_input"]')[0].attrib.get('placeholder')
    assert expected_textarea_placeholder == textarea_placeholder, 'Got %s' % textarea_placeholder


@then(r'there are radio buttons labeled "(?P<expected_sequence_type_label>[^"]*)"')
def check_sequence_type_label(context, expected_sequence_type_label):
    """
    tests the presence of radio button
    :param context: behave context
    :param expected_sequence_type_label: expected label
    """
    sequence_type_label = ''
    if TEST == 'acceptance':
        sequence_type_label = context.browser.find_elements_by_css_selector('label[for="id_seq_type_0"]')[0].text
    elif TEST == 'functional':
        sequence_type_label = context.display.cssselect('input[id="id_seq_type_0"]')[0].label.text_content()
    assert expected_sequence_type_label == sequence_type_label, 'Got %s' % sequence_type_label


@then(r'there is a "(?P<expected_radio_button_text>[^"]*)" (?P<choice>[^"]*) type radio button')
def check_sequence_type_radio_buttons(context, expected_radio_button_text, choice):
    """
    tests that the text of the radio button is as expected
    :param context: behave context
    :param expected_radio_button_text: expected radio button text
    :param choice: radio button type
    """
    if TEST == 'acceptance':
        radio_button = {
            'input sequence': {
                'DNA': context.browser.find_elements_by_css_selector('label[for="id_seq_type_1"]')[0].text,
                'Protein': context.browser.find_elements_by_css_selector('label[for="id_seq_type_0"]')[1].text
            },
            'consensus': {
                'Identity': context.browser.find_elements_by_css_selector('label[for="id_cons_type_0"]')[1].text,
                'Substitution Matrix': context.browser.find_elements_by_css_selector('label[for="id_cons_type_1"]')[
                    0].text
            },
        }
        radio_button_text = radio_button[choice][expected_radio_button_text]
        assert expected_radio_button_text == radio_button_text, 'Got %s' % radio_button_text
    elif TEST == 'functional':
        radio_button = {
            'input sequence': {
                'DNA': context.display.cssselect('label[for="id_seq_type_1"]')[0].text_content(),
                'Protein': context.display.cssselect('label[for="id_seq_type_0"]')[1].text_content()
            },
            'consensus': {
                'Identity': context.display.cssselect('label[for="id_cons_type_0"]')[1].text_content(),
                'Substitution Matrix': context.display.cssselect('label[for="id_cons_type_1"]')[0].text_content()
            },
        }
        radio_button_text = radio_button[choice][expected_radio_button_text]
        assert re.match(
            '^\s+' + re.escape(expected_radio_button_text) + '$',
            radio_button_text
        ), 'Got %s' % radio_button_text


@then(r'the "(?P<button_type>[^"]*)" button is checkable')
def check_type_buttons_are_checkable(context, button_type):
    """
    tests that the radio buttons are checkable
    :param context: behave context
    :param button_type: button to check
    """
    radio_button = {
        'DNA': {'type': 'seq', 'number': '1'},
        'Protein': {'type': 'seq', 'number': '0'},
        'Identity': {'type': 'cons', 'number': '0'},
        'Substitution Matrix': {'type': 'cons', 'number': '1'}
    }
    button_checkable = context.display.cssselect(
        'input[id="id_%s_type_%s"]' % (
            radio_button[button_type]['type'],
            radio_button[button_type]['number']
        )
    )[0].checkable
    assert button_checkable, 'Got %s' % button_checkable


@then(r'the active radio buttons are checked')
def check_active_radio_button_is_checked(context):
    """
    tests that active radio buttons are checked when clicked
    :param context: behave context
    """
    if TEST == 'acceptance':
        for radio_button in context.active_radio_buttons:
            button = context.browser.find_element_by_css_selector(
                'input[id="id_%s_type_%s"]' % (
                    radio_button['type'],
                    radio_button['number']
                )
            )
            button.click()
            button_checked = button.is_selected()
            assert button_checked, 'Expected: True\nGot %s for id_%s_type_%s' % (
                button_checked,
                radio_button['type'],
                radio_button['number']
            )
    elif TEST == 'functional':
        for radio_button in context.active_radio_buttons:
            button_checked = context.display.cssselect(
                'input[id="id_%s_type_%s"]' % (
                    radio_button['type'],
                    radio_button['number']
                )
            )[0].checked
            assert button_checked, 'Expected: True\nGot %s for id_%s_type_%s' % (
                button_checked,
                radio_button['type'],
                radio_button['number']
            )


@then(r'the inactive radio buttons are not checked')
def check_active_radio_button_is_checked(context):
    """
    tests that inactive radio buttons are not checked when active buttons are clicked
    :param context: behave context
    """
    if TEST == 'acceptance':
        for i, radio_button in enumerate(context.active_radio_buttons):
            button = context.browser.find_element_by_css_selector(
                'input[id="id_%s_type_%s"]' % (
                    radio_button['type'],
                    radio_button['number']
                )
            )
            button.click()
            button_checked = context.browser.find_element_by_css_selector(
                'input[id="id_%s_type_%s"]' % (
                    context.inactive_radio_buttons[i]['type'],
                    context.inactive_radio_buttons[i]['number']
                )
            ).is_selected()
            assert not button_checked, 'Expected: False\nGot %s for id_%s_type_%s' % (
                button_checked,
                radio_button['type'],
                radio_button['number']
            )


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
    if TEST == 'acceptance':
        button = context.browser.find_element_by_css_selector('button[id="%s"]' % button_id[button_text])
        assert button_text == button.text, 'Got %s' % button.text
        assert button.get_attribute('name') == button_name, 'Got %s' % button.get_attribute('name')
        assert button.get_attribute('value') == button_value, 'Got %s' % button.get_attribute('value')
    elif TEST == 'functional':
        button = context.display.cssselect('button[id="%s"]' % button_id[button_text])
        assert re.match(
            '^\\n\s+' + re.escape(button_text) + '\\n\s+$', button[0].text_content()
        ), 'Got %s' % button[0].text_content()
        assert button[0].attrib.get('name') == button_name, 'Got %s' % button[0].attrib.get('name')
        assert button[0].attrib.get('value') == button_value, 'Got %s' % button[0].attrib.get('value')


@then(r'there is a "(?P<link_text>[^"]*)" button with "(?P<href>[^"]*)" href')
def check_link_href(context, link_text, href):
    """
    tests the presence and attributes of a link
    :param context: behave context
    :param link_text: link text
    :param href: URL of link
    """
    if TEST == 'acceptance':
        home_button = context.browser.find_elements_by_css_selector('a[class="navbar-brand"]')

        assert link_text == home_button[0].text, 'Expected: %s\nGot: %s' % (link_text, home_button[0].text)
        browser_href = ''
        if context.browser_name == 'Firefox':
            if SERVER_URL == 'liveserver':
                href_expected = href
            else:
                href_expected = SERVER_URL + href
            browser_href = href_expected
        elif context.browser_name == 'Chrome':
            if SERVER_URL == 'liveserver':
                href_expected = context.base_url + href
            else:
                href_expected = SERVER_URL + href
            browser_href = href_expected
        assert browser_href == home_button[0].get_attribute('href'), \
            '\nExpected: %s\nGot: %s' % (browser_href, home_button[0].get_attribute('href'))
    elif TEST == 'functional':
        home_button = context.display.cssselect('a[class="navbar-brand"]')
        assert link_text == home_button[0].text_content(), 'Got %s' % home_button[0].text_content()
        assert href == home_button[0].attrib.get('href'), \
            '\nExpected: %s\nGot: %s' % (href, home_button[0].attrib.get('href'))
