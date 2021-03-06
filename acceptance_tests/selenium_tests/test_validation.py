import time
from formalign.settings import CHROME_DRIVER, SERVER_URL, BROWSER_DELAY, TEST_CASE, FIREFOX_BINARY
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.desired_capabilities import DesiredCapabilities
import pyperclip
from helper_funcs.helpers_test import file_to_string
from base.forms import EMPTY_ERROR, FORMAT_ERROR, CHARACTER_ERROR, ALIGNMENT_ERROR, LESS_THAN_TWO_SEQS_ERROR


class InputValidationTestCaseChrome(TEST_CASE):
    def setUp(self):
        self.browser = webdriver.Chrome(
            CHROME_DRIVER
        )
        self.delay = float(BROWSER_DELAY)
        if SERVER_URL == 'liveserver':
            self.url = self.live_server_url
            self.wait = self.delay * 5
        else:
            self.url = SERVER_URL
            self.wait = self.delay
        self.browser.implicitly_wait(5)

    def tearDown(self):
        self.browser.quit()

    def invalid_format_sequence(self, file):
        seq_type_button_dict = {'DNA': 'input#id_seq_type_1', 'protein': 'input#id_seq_type_0'}
        test_seq = {'DNA': 'AGTCC-TAAGGTCGCCAATGGGCA', 'protein': 'MKERBGWAQ--QGKKPWRF--EEW'}

        # she visits the Formalign.eu site
        self.browser.get(self.url + '/')

        # She clicks the appropriate button and clears the input field
        seq_type_button = self.browser.find_element_by_css_selector(seq_type_button_dict[file['seq_type']])
        seq_type_button.click()
        self.assertEqual(
            True,
            seq_type_button.is_selected(),
            'button is selected for ' + file['align_format'] + ' ' + file['seq_type'] + ': ' +
            str(seq_type_button.is_selected())
        )
        alignment_input = self.browser.find_element_by_css_selector('textarea#id_align_input')
        alignment_input.clear()

        # She pastes in an invalid fasta alignment
        alignment_string = file_to_string(file['invalid'])
        pyperclip.copy(alignment_string)
        alignment_input.send_keys(Keys.CONTROL, 'v')
        self.browser.find_element_by_id('submit-align').click()

        # Wait for Firefox
        time.sleep(self.wait)

        # Since her FASTA format is invalid she gets redirected to the submission form where she sees an
        # error message telling her that her alignment format is invalid
        self.assertEqual(
            'Formalign.eu Home',
            self.browser.title,
            'browser.title for ' + file['align_format'] + ' ' + file['seq_type'] + ': ' + self.browser.title
        )
        error = self.browser.find_element_by_css_selector('.errorlist').find_element_by_tag_name('li')
        self.assertEqual(
            FORMAT_ERROR,
            error.text,
            'error.text for ' + file['align_format'] + ' ' + file['seq_type'] + ': ' + error.text,

        )

        # she corrects her alignment and submits an invalid clustal alignment
        alignment_input = self.browser.find_element_by_css_selector('textarea#id_align_input')
        alignment_input.clear()
        alignment_string = file_to_string(file['valid'])
        pyperclip.copy(alignment_string)
        alignment_input.send_keys(Keys.CONTROL, 'v')
        self.browser.find_element_by_id('submit-align').click()

        # Wait for Firefox
        time.sleep(self.wait * 10)

        # She got it right this time and is redirected to a page showing the submitted sequences from her alignment
        self.assertEqual(
            'Formalign.eu Sequence Display',
            self.browser.title,
            'browser.title for ' + file['align_format'] + ' ' + file['seq_type'] + ': ' + self.browser.title)
        first_seq_info = self.browser.find_elements_by_css_selector('.query_seq_meta')[0]
        self.assertEqual(
            'sequence1:',
            first_seq_info.text,
            'seq id for ' + file['align_format'] + ' ' + file['seq_type'] + ': ' + first_seq_info.text

        )
        first_seq_content = self.browser.find_elements_by_css_selector('.query_seq_display')[0]
        self.assertIsNotNone(first_seq_content)
        self.assertEqual(
            test_seq[file['seq_type']],
            first_seq_content.text,
            'seq id for ' + file['align_format'] + ' ' + file['seq_type'] + ': ' + first_seq_content.text
        )

        # She wonders whether she can use other formats and decides to navigate back to the home page
        home_button = self.browser.find_element_by_css_selector('.navbar-brand')
        home_button.click()

    def alignment_validation(self, **kwargs):
        seq_type_button_dict = {'DNA': 'input#id_seq_type_1', 'protein': 'input#id_seq_type_0'}
        test_seq = {'DNA': 'AGTCC-TAAGGTCGCCAATGGGCA', 'protein': 'MKERBGWAQ--QGKKPWRF--EEW'}

        # User visits the formalign.eu site.
        self.browser.get(self.url + '/')

        # She clicks the DNA button
        seq_type_button = self.browser.find_element_by_css_selector(seq_type_button_dict[kwargs['seq_type']])
        seq_type_button.click()

        # She is so excited that she inadvertently submits the empty input
        self.browser.find_element_by_id('submit-align').click()

        # Wait for Firefox
        time.sleep(self.wait)

        self.assertEqual('Formalign.eu Home', self.browser.title, self.browser.title)
        error = self.browser.find_element_by_css_selector('.errorlist').find_element_by_tag_name('li')
        self.assertEqual(
            EMPTY_ERROR,
            error.text
        )

        # She decides to try it out so she pastes in an alignment and submits
        alignment_input = self.browser.find_element_by_css_selector('textarea#id_align_input')
        alignment_string = file_to_string(kwargs['seq_type'] + '_invalid_fasta.fasta')
        alignment_input.send_keys(alignment_string)
        self.browser.find_element_by_id('submit-align').click()

        # Wait for Firefox
        time.sleep(self.wait * 5)

        # unfortunately her FASTA format is invalid so she gets redirected to the submission form where she sees an
        # error message telling her that her FASTA format is invalid
        self.assertEqual(self.browser.title, 'Formalign.eu Home', self.browser.title)
        error = self.browser.find_element_by_css_selector('.errorlist').find_element_by_tag_name('li')
        self.assertEqual(
            FORMAT_ERROR,
            error.text
        )

        # she corrects her alignment and resubmits
        alignment_string = file_to_string(kwargs['seq_type'] + '_invalid_characters.fasta')
        alignment_input = self.browser.find_element_by_css_selector('textarea#id_align_input')
        alignment_input.clear()
        alignment_input.send_keys(alignment_string)
        self.browser.find_element_by_id('submit-align').click()

        # Wait for Firefox
        time.sleep(self.wait * 2)

        # unfortunately now her sequences contain invalid characters so she gets redirected to the submission form
        # again where she sees an error message telling her that her sequences contain invalid characters
        self.assertEqual(self.browser.title, 'Formalign.eu Home', self.browser.title)
        error = self.browser.find_element_by_css_selector('.errorlist').find_element_by_tag_name('li')
        self.assertEqual(
            CHARACTER_ERROR + 'sequence1',
            error.text

        )

        # she corrects her alignment again and resubmits
        alignment_string = file_to_string(kwargs['seq_type'] + '_too_few_sequences.fasta')
        alignment_input = self.browser.find_element_by_css_selector('textarea#id_align_input')
        alignment_input.clear()
        alignment_input.send_keys(alignment_string)
        self.browser.find_element_by_id('submit-align').click()

        # Wait for Firefox
        time.sleep(self.wait * 2)

        # unfortunately this time she accidentally erased one sequence and is left with only one sequence so she gets
        # redirected to the submission form again where she sees an error message telling her that her alignment is not
        # an alignment since it contains only one sequence
        self.assertEqual('Formalign.eu Home', self.browser.title, self.browser.title)
        error = self.browser.find_element_by_css_selector('.errorlist').find_element_by_tag_name('li')
        self.assertEqual(
            LESS_THAN_TWO_SEQS_ERROR,
            error.text,
        )

        # she adds the missing sequence and resubmits
        alignment_string = file_to_string(kwargs['seq_type'] + '_invalid_alignment.fasta')
        alignment_input = self.browser.find_element_by_css_selector('textarea#id_align_input')
        alignment_input.clear()
        alignment_input.send_keys(alignment_string)
        self.browser.find_element_by_id('submit-align').click()

        # Wait for Firefox
        time.sleep(self.wait)

        # it must be starting to be a bit late since she added some residues to her first sequence so it is longer than
        # the second now so she gets redirected to the submission form again where she sees an error message telling her
        # that her alignment is not an alignment since the sequences do not all have the same length
        self.assertEqual(self.browser.title, 'Formalign.eu Home', self.browser.title)
        error = self.browser.find_element_by_css_selector('.errorlist').find_element_by_tag_name('li')
        self.assertEqual(
            ALIGNMENT_ERROR,
            error.text,
        )

        # She tries one final time and threatens to throw her laptop out of the window if she gets another
        # error message
        alignment_string = file_to_string(kwargs['seq_type'] + '.fasta')
        alignment_input = self.browser.find_element_by_css_selector('textarea#id_align_input')
        alignment_input.clear()
        alignment_input.send_keys(alignment_string)
        self.browser.find_element_by_id('submit-align').click()

        # Wait for Firefox
        time.sleep(self.wait * 5)

        # She got it right this time and is redirected to a page showing the submitted sequences from her alignment
        self.assertEqual(self.browser.title, 'Formalign.eu Sequence Display', self.browser.title)
        first_seq_info = self.browser.find_elements_by_css_selector('.query_seq_meta')[0]
        self.assertEqual(
            'sequence1:',
            first_seq_info.text
        )
        first_seq_content = self.browser.find_elements_by_css_selector('.query_seq_display')[0]
        self.assertIsNotNone(first_seq_content)
        self.assertEqual(first_seq_content.text, test_seq[kwargs['seq_type']])

    def test_alignment_validation(self):
        tests = [{'seq_type': 'DNA'}, {'seq_type': 'protein'}]
        for t in tests:
            self.alignment_validation(**t)

    def test_alignment_format_validation(self):
        # User visits the formalign.eu site.
        # self.browser.get(self.url + '/')

        # She tries a number of invalid and valid alignment formats
        files = [
            {'valid': 'protein.fasta', 'invalid': 'protein_invalid_fasta.fasta', 'seq_type': 'protein',
             'align_format': 'fasta'},
            {'valid': 'DNA.fasta', 'invalid': 'DNA_invalid_fasta.fasta', 'seq_type': 'DNA',
             'align_format': 'fasta'},
        ]
        for file in files:
            self.invalid_format_sequence(file)


class InputValidationTestCaseFirefox(InputValidationTestCaseChrome):
    def setUp(self):
        caps = DesiredCapabilities.FIREFOX
        caps['marionette'] = True
        caps['binary'] = FIREFOX_BINARY
        self.browser = webdriver.Firefox(capabilities=caps)
        self.delay = float(BROWSER_DELAY)
        if SERVER_URL == 'liveserver':
            self.url = self.live_server_url
            self.wait = self.delay * 10
        else:
            self.url = SERVER_URL
            self.wait = self.delay * 2.5

    def tearDown(self):
        self.browser.quit()
