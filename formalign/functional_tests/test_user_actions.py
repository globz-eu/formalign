from django.contrib.staticfiles.testing import StaticLiveServerTestCase
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
import pyperclip
import time
from helper_funcs.helpers_test import file_to_string

__author__ = 'Stefan Dieterle'


class BasicUserTestCase(StaticLiveServerTestCase):
    def setUp(self):
        # self.browser = webdriver.Chrome(
        #     '/home/golgotux/Dropbox/Documents/Backups/Python Scripts/pycharm_helpers/chromedriver'
        # )
        self.browser = webdriver.Firefox()
        self.browser.implicitly_wait(2)

    def tearDown(self):
        # time.sleep(5)
        self.browser.quit()

    def test_basic_user_experience(self):
        """
        Tests basic user interaction with formalign.eu site
        :return:
        """
        # Lambda user is a biologist who has to make a nice figure containing a multiple alignment for a presentation.
        # She visits the formalign.eu site.
        self.browser.get(self.live_server_url + '/')

        # User sees she's on the right page because she can see the name of the site in the heading.
        self.assertEqual(self.browser.title, 'Formalign.eu Home', self.browser.title)
        brand_element = self.browser.find_element_by_css_selector('.navbar-brand')
        self.assertEqual('Formalign.eu', brand_element.text)

        # She sees a form that says 'Paste in your alignment in FASTA format:'
        alignment_input = self.browser.find_element_by_css_selector('textarea#id_align_input')
        self.assertIsNotNone(self.browser.find_element_by_css_selector('label[for="id_align_input"]'))
        self.assertEqual(alignment_input.get_attribute('placeholder'), 'FASTA alignment')

        # She sees two radio buttons for DNA and protein
        dna_button = self.browser.find_element_by_css_selector('input#id_seq_type_1')
        self.assertIsNotNone(dna_button)
        protein_button = self.browser.find_element_by_css_selector('input#id_seq_type_0')
        self.assertIsNotNone(protein_button)

        # She sees that the DNA button is selected by default
        self.assertEqual(dna_button.is_selected(), True)

        # She clicks the Protein radio button and sees that it gets selected and the DNA button gets unselected
        protein_button.click()
        self.assertEqual(protein_button.is_selected(), True)
        self.assertEqual(dna_button.is_selected(), False)

        # She pastes in a protein alignment to see what happens
        alignment_string = file_to_string('spa_align_clustal_omega.fasta')
        pyperclip.copy(alignment_string)
        alignment_input = self.browser.find_element_by_css_selector('textarea#id_align_input')
        alignment_input.send_keys(Keys.CONTROL, 'v')
        self.browser.find_element_by_id('submit-fasta').click()

        # She is redirected to a page showing the submitted sequences from her alignment
        self.assertEqual(self.browser.title, 'Formalign.eu Sequence Display', self.browser.title)
        seq_content = self.browser.find_elements_by_css_selector('.query_seq_display')
        self.assertIsNotNone(seq_content)
        for f in seq_content:
            self.assertTrue(len(f.text) <= 80)
        self.fail('Incomplete Test')

    def test_DNA_alignment_user_experience(self):
        # User visits the formalign.eu site.
        self.browser.get(self.live_server_url + '/')

        # She clicks the DNA button
        dna_button = self.browser.find_element_by_css_selector('input#id_seq_type_1')
        dna_button.click()

        # She decides to try the default with a DNA alignment first so she pastes in a DNA alignment and submits
        alignment_input = self.browser.find_element_by_css_selector('textarea#id_align_input')
        alignment_string = file_to_string('DNA.fasta')
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

    def test_protein_alignment_user_experience(self):
        # User visits the formalign.eu site.
        self.browser.get(self.live_server_url + '/')

        # she clicks the protein radio button
        protein_button = self.browser.find_element_by_css_selector('input#id_seq_type_0')
        protein_button.click()

        # She types in an alignment of her favorite proteins and submits it.
        alignment_input = self.browser.find_element_by_css_selector('textarea#id_align_input')
        alignment_string = file_to_string('short_invalid_fasta.fasta')
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
        alignment_string = file_to_string('short_invalid_characters.fasta')
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
        alignment_string = file_to_string('short_too_few_sequences.fasta')
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
        alignment_string = file_to_string('short_invalid_alignment.fasta')
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
        alignment_string = file_to_string('short.fasta')
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
