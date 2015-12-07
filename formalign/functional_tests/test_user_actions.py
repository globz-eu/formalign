from django.contrib.staticfiles.testing import StaticLiveServerTestCase
from selenium import webdriver
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
        # time.sleep(2)
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

        # She decides to give it a try. She types in an alignment of her favorite proteins and submits it.
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
