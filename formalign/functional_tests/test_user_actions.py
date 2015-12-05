from django.contrib.staticfiles.testing import StaticLiveServerTestCase
from formalign.settings import BASE_DIR
from selenium import webdriver
import time
import os

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
        home_page = self.browser.get(self.live_server_url + '/')

        # User sees she's on the right page because she can see the name of the site in the heading.
        self.assertEqual(self.browser.title, 'Formalign.eu Home', self.browser.title)
        brand_element = self.browser.find_element_by_css_selector('.navbar-brand')
        self.assertEqual('Formalign.eu', brand_element.text)

        # She sees a form that says 'Paste in your alignment in FASTA format:'
        alignment_input = self.browser.find_element_by_css_selector('textarea#id_align_input')
        self.assertIsNotNone(self.browser.find_element_by_css_selector('label[for="id_align_input"]'))
        self.assertEqual(alignment_input.get_attribute('placeholder'), 'FASTA alignment')

        # She decides to give it a try. She pastes in her alignment of favorite proteins and submits it.
        with open(os.path.join(BASE_DIR, 'test_data/short.fasta'), 'r') as input_seqs:
            alignment_string = input_seqs.read()
        alignment_input.send_keys(alignment_string)
        self.browser.find_element_by_id('submit-fasta').click()

        # She is redirected to a page showing the submitted sequences from her alignment
        self.assertEqual(self.browser.title, 'Formalign.eu Sequence Display', self.browser.title)
        first_seq_info = self.browser.find_elements_by_css_selector('.query_seq_meta')[0]
        self.assertEqual(
            first_seq_info.text,
            'Short sequence1'
        )
        first_seq_content = self.browser.find_elements_by_css_selector('.query_seq_display')[0]
        self.assertIsNotNone(first_seq_content)
        self.assertEqual(first_seq_content.text, 'MKERBGWAQ--QGKKPWRF--EEW')
        self.fail('Incomplete Test')

        # She is redirected to a display page where she sees her alignment rendered in the default way.
