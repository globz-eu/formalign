from django.contrib.staticfiles.testing import StaticLiveServerTestCase
from private import CHROME_DRIVER
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
import pyperclip
from helper_funcs.helpers_test import file_to_string

__author__ = 'Stefan Dieterle'


class InputValidationTestCase(StaticLiveServerTestCase):
    def setUp(self):
        # self.browser = webdriver.Chrome(
        #     CHROME_DRIVER
        # )
        self.browser = webdriver.Firefox()
        self.browser.implicitly_wait(2)

    def tearDown(self):
        # time.sleep(5)
        self.browser.quit()

    def invalid_format_test_sequence(self, **kwargs):
        seq_type_button = {'DNA': 'input#id_seq_type_1', 'protein': 'input#id_seq_type_0'}

        # She clicks the appropriate button and clears the input field
        protein_button = self.browser.find_element_by_css_selector(seq_type_button[kwargs['seq_type']])
        protein_button.click()
        self.assertEqual(
                True,
                protein_button.is_selected(),
                'button is selected for ' + kwargs['align_format'] + ' ' + kwargs['seq_type'] + ': ' +
                str(protein_button.is_selected())
        )
        alignment_input = self.browser.find_element_by_css_selector('textarea#id_align_input')
        alignment_input.clear()

        # She pastes in an invalid fasta alignment
        alignment_string = file_to_string(kwargs['invalid'])
        pyperclip.copy(alignment_string)
        alignment_input.send_keys(Keys.CONTROL, 'v')
        self.browser.find_element_by_id('submit-align').click()

        # Since her FASTA format is invalid she gets redirected to the submission form where she sees an
        # error message telling her that her alignment format is invalid
        self.assertEqual(
                'Formalign.eu Home',
                self.browser.title,
                'browser.title for ' + kwargs['align_format'] + ' ' + kwargs['seq_type'] + ': ' + self.browser.title
        )
        error = self.browser.find_element_by_css_selector('.errorlist').find_element_by_tag_name('li')
        self.assertEqual(
                'The server could not figure out what format this is, '
                'please double check your input or try a different format',
                error.text,
                'error.text for ' + kwargs['align_format'] + ' ' + kwargs['seq_type'] + ': ' + error.text,

        )

        # she corrects her alignment and submits an invalid clustal alignment
        alignment_input = self.browser.find_element_by_css_selector('textarea#id_align_input')
        alignment_input.clear()
        alignment_string = file_to_string(kwargs['valid'])
        pyperclip.copy(alignment_string)
        alignment_input.send_keys(Keys.CONTROL, 'v')
        self.browser.find_element_by_id('submit-fasta').click()

        # She got it right this time and is redirected to a page showing the submitted sequences from her alignment
        self.assertEqual(
                'Formalign.eu Sequence Display',
                self.browser.title,
                'browser.title for ' + kwargs['align_format'] + ' ' + kwargs['seq_type'] + ': ' + self.browser.title)
        first_seq_info = self.browser.find_elements_by_css_selector('.query_seq_meta')[0]
        self.assertEqual(
                'sequence1:',
                first_seq_info.text,
                'seq id for ' + kwargs['align_format'] + ' ' + kwargs['seq_type'] + ': ' + first_seq_info.text

        )
        first_seq_content = self.browser.find_elements_by_css_selector('.query_seq_display')[0]
        self.assertIsNotNone(first_seq_content)
        self.assertEqual(
                'MKERBGWAQ--QGKKPWRF--EEW',
                first_seq_content.text,
                'seq id for ' + kwargs['align_format'] + ' ' + kwargs['seq_type'] + ': ' + first_seq_content.text
        )

        # She wonders whether she can use other formats and decides to navigate back to the home page
        home_button = self.browser.find_element_by_css_selector('.navbar-brand')
        home_button.click()

    def test_alignment_format_validation(self):
        # User visits the formalign.eu site.
        self.browser.get(self.live_server_url + '/')

        # She tries a number of invalid and valid alignment formats
        files = [
            {'valid': 'protein.fasta', 'invalid': 'protein_invalid_fasta.fasta', 'seq_type': 'protein',
             'align_format': 'fasta'},
            {'valid': 'protein.clustal', 'invalid': 'protein_invalid_clustal.clustal', 'seq_type': 'protein',
             'align_format': 'clustal'},
            {'valid': 'protein.ig', 'invalid': 'protein_invalid_ig.ig', 'seq_type': 'protein',
             'align_format': 'ig'},
            {'valid': 'protein.nexus', 'invalid': 'protein_invalid_nexus.nexus', 'seq_type': 'protein',
             'align_format': 'nexus'},
            {'valid': 'protein.phylip', 'invalid': 'protein_invalid_phylip.phylip', 'seq_type': 'protein',
             'align_format': 'phylip'},
            {'valid': 'protein.stockholm', 'invalid': 'protein_invalid_stockholm.sto', 'seq_type': 'protein',
             'align_format': 'stockholm'},
            {'valid': 'DNA.fasta', 'invalid': 'DNA_invalid_fasta.fasta', 'seq_type': 'DNA',
             'align_format': 'fasta'},
            {'valid': 'DNA.clustal', 'invalid': 'DNA_invalid_clustal.clustal', 'seq_type': 'DNA',
             'align_format': 'clustal'},
            {'valid': 'DNA.ig', 'invalid': 'DNA_invalid_ig.ig', 'seq_type': 'DNA',
             'align_format': 'ig'},
            {'valid': 'DNA.nexus', 'invalid': 'DNA_invalid_nexus.nexus', 'seq_type': 'DNA',
             'align_format': 'nexus'},
            {'valid': 'DNA.phylip', 'invalid': 'DNA_invalid_phylip.phylip', 'seq_type': 'DNA',
             'align_format': 'phylip'},
            {'valid': 'DNA.stockholm', 'invalid': 'DNA_invalid_stockholm.sto', 'seq_type': 'DNA',
             'align_format': 'stockholm'},

        ]
        for file in files:
            self.invalid_format_test_sequence(**file)
        self.fail('Incomplete Test')

    def test_DNA_alignment_validation(self):
        # User visits the formalign.eu site.
        self.browser.get(self.live_server_url + '/')

        # She clicks the DNA button
        dna_button = self.browser.find_element_by_css_selector('input#id_seq_type_1')
        dna_button.click()

        # She decides to try the default with a DNA alignment so she pastes in a DNA alignment and submits
        alignment_input = self.browser.find_element_by_css_selector('textarea#id_align_input')
        alignment_string = file_to_string('DNA_invalid_fasta.fasta')
        alignment_input.send_keys(alignment_string)
        self.browser.find_element_by_id('submit-fasta').click()

        # unfortunately her FASTA format is invalid so she gets redirected to the submission form where she sees an
        # error message telling her that her FASTA format is invalid
        self.assertEqual(self.browser.title, 'Formalign.eu Home', self.browser.title)
        error = self.browser.find_element_by_css_selector('.errorlist').find_element_by_tag_name('li')
        self.assertEqual(
                error.text,
                'Sequence is not FASTA compliant, no ">" as first character'
        )

        # she corrects her alignment and resubmits
        alignment_string = file_to_string('DNA_invalid_characters.fasta')
        alignment_input = self.browser.find_element_by_css_selector('textarea#id_align_input')
        alignment_input.clear()
        alignment_input.send_keys(alignment_string)
        self.browser.find_element_by_id('submit-fasta').click()

        # unfortunately now her sequences contain invalid characters so she gets redirected to the submission form
        # again where she sees an error message telling her that her sequences contain invalid characters
        self.assertEqual(self.browser.title, 'Formalign.eu Home', self.browser.title)
        error = self.browser.find_element_by_css_selector('.errorlist').find_element_by_tag_name('li')
        self.assertEqual(
                error.text,
                'Invalid character in sequence: sequence1'
        )

        # she corrects her alignment again and resubmits
        alignment_string = file_to_string('DNA_too_few_sequences.fasta')
        alignment_input = self.browser.find_element_by_css_selector('textarea#id_align_input')
        alignment_input.clear()
        alignment_input.send_keys(alignment_string)
        self.browser.find_element_by_id('submit-fasta').click()

        # unfortunately this time she accidentally erased one sequence and is left with only one sequence so she gets
        # redirected to the submission form again where she sees an error message telling her that her alignment is not
        # an alignment since it contains only one sequence
        self.assertEqual(self.browser.title, 'Formalign.eu Home', self.browser.title)
        error = self.browser.find_element_by_css_selector('.errorlist').find_element_by_tag_name('li')
        self.assertEqual(
                error.text,
                'Submitted data is not a valid alignment, it contains less than 2 sequences'
        )

        # she adds the missing sequence and resubmits
        alignment_string = file_to_string('DNA_invalid_alignment.fasta')
        alignment_input = self.browser.find_element_by_css_selector('textarea#id_align_input')
        alignment_input.clear()
        alignment_input.send_keys(alignment_string)
        self.browser.find_element_by_id('submit-fasta').click()

        # it must be starting to be a bit late since she added some residues to her first sequence so it is longer than
        # the second now so she gets redirected to the submission form again where she sees an error message telling her
        # that her alignment is not an alignment since the sequences do not all have the same length
        self.assertEqual(self.browser.title, 'Formalign.eu Home', self.browser.title)
        error = self.browser.find_element_by_css_selector('.errorlist').find_element_by_tag_name('li')
        self.assertEqual(
                error.text,
                'Alignment invalid, sequences have different lengths'
        )

        # She tries one final time and threatens to throw her laptop out of the window if she gets another
        # error message
        alignment_string = file_to_string('DNA.fasta')
        alignment_input = self.browser.find_element_by_css_selector('textarea#id_align_input')
        alignment_input.clear()
        alignment_input.send_keys(alignment_string)
        self.browser.find_element_by_id('submit-fasta').click()

        # She got it right this time and is redirected to a page showing the submitted sequences from her alignment
        self.assertEqual(self.browser.title, 'Formalign.eu Sequence Display', self.browser.title)
        first_seq_info = self.browser.find_elements_by_css_selector('.query_seq_meta')[0]
        self.assertEqual(
                first_seq_info.text,
                'sequence1:'
        )
        first_seq_content = self.browser.find_elements_by_css_selector('.query_seq_display')[0]
        self.assertIsNotNone(first_seq_content)
        self.assertEqual(first_seq_content.text, 'AGTCC-TAAGGTCGCCAATGGGCA')

        # She is redirected to a display page where she sees her alignment rendered in the default way.
        self.fail('Incomplete Test')

    def test_protein_alignment_validation(self):
        # User visits the formalign.eu site.
        self.browser.get(self.live_server_url + '/')

        # she clicks the protein radio button
        protein_button = self.browser.find_element_by_css_selector('input#id_seq_type_0')
        protein_button.click()

        # She types in an alignment of her favorite proteins and submits it.
        alignment_input = self.browser.find_element_by_css_selector('textarea#id_align_input')
        alignment_string = file_to_string('protein_invalid_fasta.fasta')
        alignment_input.send_keys(alignment_string)
        self.browser.find_element_by_id('submit-fasta').click()

        # unfortunately her FASTA format is invalid so she gets redirected to the submission form where she sees an
        # error message telling her that her FASTA format is invalid
        self.assertEqual(self.browser.title, 'Formalign.eu Home', self.browser.title)
        error = self.browser.find_element_by_css_selector('.errorlist').find_element_by_tag_name('li')
        self.assertEqual(
                error.text,
                'Sequence is not FASTA compliant, no ">" as first character'
        )

        # she corrects her alignment and resubmits
        alignment_string = file_to_string('protein_invalid_characters.fasta')
        alignment_input = self.browser.find_element_by_css_selector('textarea#id_align_input')
        alignment_input.clear()
        alignment_input.send_keys(alignment_string)
        self.browser.find_element_by_id('submit-fasta').click()

        # unfortunately now her sequences contain invalid characters so she gets redirected to the submission form
        # again where she sees an error message telling her that her sequences contain invalid characters
        self.assertEqual(self.browser.title, 'Formalign.eu Home', self.browser.title)
        error = self.browser.find_element_by_css_selector('.errorlist').find_element_by_tag_name('li')
        self.assertEqual(
                error.text,
                'Invalid character in sequence: Short sequence3'
        )

        # she corrects her alignment again and resubmits
        alignment_string = file_to_string('protein_too_few_sequences.fasta')
        alignment_input = self.browser.find_element_by_css_selector('textarea#id_align_input')
        alignment_input.clear()
        alignment_input.send_keys(alignment_string)
        self.browser.find_element_by_id('submit-fasta').click()

        # unfortunately this time she accidentally erased one sequence and is left with only one sequence so she gets
        # redirected to the submission form again where she sees an error message telling her that her alignment is not
        # an alignment since it contains only one sequence
        self.assertEqual(self.browser.title, 'Formalign.eu Home', self.browser.title)
        error = self.browser.find_element_by_css_selector('.errorlist').find_element_by_tag_name('li')
        self.assertEqual(
                error.text,
                'Submitted data is not a valid alignment, it contains less than 2 sequences'
        )

        # she adds the missing sequence and resubmits
        alignment_string = file_to_string('protein_invalid_alignment.fasta')
        alignment_input = self.browser.find_element_by_css_selector('textarea#id_align_input')
        alignment_input.clear()
        alignment_input.send_keys(alignment_string)
        self.browser.find_element_by_id('submit-fasta').click()

        # it must be starting to be a bit late since she added some residues to her first sequence so it is longer than
        # the second now so she gets redirected to the submission form again where she sees an error message telling her
        # that her alignment is not an alignment since the sequences do not all have the same length
        self.assertEqual(self.browser.title, 'Formalign.eu Home', self.browser.title)
        error = self.browser.find_element_by_css_selector('.errorlist').find_element_by_tag_name('li')
        self.assertEqual(
                error.text,
                'Alignment invalid, sequences have different lengths'
        )

        # She tries one final time and threatens to throw her laptop out of the window if she gets another
        # error message
        alignment_string = file_to_string('protein.fasta')
        alignment_input = self.browser.find_element_by_css_selector('textarea#id_align_input')
        alignment_input.clear()
        alignment_input.send_keys(alignment_string)
        self.browser.find_element_by_id('submit-fasta').click()

        # She got it right this time and is redirected to a page showing the submitted sequences from her alignment
        self.assertEqual(self.browser.title, 'Formalign.eu Sequence Display', self.browser.title)
        first_seq_info = self.browser.find_elements_by_css_selector('.query_seq_meta')[0]
        self.assertEqual(
                first_seq_info.text,
                'Short sequence1:'
        )
        first_seq_content = self.browser.find_elements_by_css_selector('.query_seq_display')[0]
        self.assertIsNotNone(first_seq_content)
        self.assertEqual(first_seq_content.text, 'MKERBGWAQ--QGKKPWRF--EEW')

        # She is redirected to a display page where she sees her alignment rendered in the default way.
        self.fail('Incomplete Test')
