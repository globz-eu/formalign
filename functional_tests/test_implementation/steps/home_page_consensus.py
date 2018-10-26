from behave import then, use_step_matcher
from formalign.settings import TEST


use_step_matcher('re')


@then(r'there is a "(?P<dropdown>[^"]*)" dropdown menu')
def dropdown_present(context, dropdown):
    dropdown_menu_label = ''
    if TEST == 'acceptance':
        context.dropdown_menu = {
            'consensus': context.browser.get_element_by_css_selector('#id_consensus_select')
        }
        dropdown_menu_label = context.dropdown_menu[dropdown].label.text
    elif TEST == 'functional':
        context.dropdown_menu = {
            'consensus': context.display.cssselect('input[id="id_consensus_select"]')[0]
        }
        dropdown_menu_label = context.dropdown_menu[dropdown].label.text_content()
    assert 'consensus' == dropdown_menu_label, 'Got %s' % dropdown_menu_label


@then(r'the "(?P<dropdown>[^"]*)" dropdown menu contains (?P<consensus_choice>.*)')
def check_dropdown_content(context, dropdown, consensus_choice):
    dropdown_menu = {'consensus': context.display.cssselect('input[id="id_consensus_select"]')[0]}
    assert consensus_choice in dropdown_menu[dropdown], 'Got %s' % dropdown_menu[dropdown]


@then(r'the "(?P<dropdown>[^"]*)" dropdown menu contents (?P<consensus_choice>.*) are hidden')
def check_dropdown_hidden_content(context, dropdown, consensus_choice):
    dropdown_menu = {'consensus': context.display.cssselect('input[id="id_consensus_select"]')[0]}
    assert consensus_choice in dropdown_menu[dropdown], 'Got %s' % dropdown_menu[dropdown]
