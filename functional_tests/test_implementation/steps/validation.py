from behave import then, use_step_matcher
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as ec
from base.forms import EMPTY_ERROR, FORMAT_ERROR, CHARACTER_ERROR, ALIGNMENT_ERROR, LESS_THAN_TWO_SEQS_ERROR
from formalign.settings import TEST


use_step_matcher('re')


@then(r'the user should see the error message: (?P<error_message>.*)')
def check_error_message(context, error_message):
    errors = {
        'empty error': EMPTY_ERROR,
        'character error': '%ssequence1' % CHARACTER_ERROR,
        'less than two sequences error': LESS_THAN_TWO_SEQS_ERROR,
        'alignment error': ALIGNMENT_ERROR,
        'format error': FORMAT_ERROR
    }
    error_text = ''
    if TEST == 'acceptance':
        try:
            WebDriverWait(context.browser, 2).until(
                ec.presence_of_element_located((By.CLASS_NAME, 'errorlist'))
            )
        finally:
            error_text = context.browser.find_elements_by_class_name(
                'errorlist'
            )[0].find_elements_by_css_selector('li')[0].text
    elif TEST == 'functional':
        error_text = context.display.cssselect('ul[class="errorlist"]')[0].cssselect('li')[0].text_content()
    assert errors[error_message] == error_text, 'Got %s' % error_text
