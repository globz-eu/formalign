from django.test import LiveServerTestCase
from formalign.settings import BASE_DIR
from selenium import webdriver
import time
import os

__author__ = 'Stefan Dieterle'


class BasicUserTestCase(LiveServerTestCase):
    def setUp(self):
        # self.browser = webdriver.Chrome(
        #     '/home/golgotux/Dropbox/Documents/Backups/Python Scripts/pycharm_helpers/chromedriver'
        # )
        self.browser = webdriver.Firefox()
        self.browser.implicitly_wait(2)

    def tearDown(self):
        time.sleep(2)
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
        brand_element = self.browser.find_element_by_css_selector('.navbar-brand')
        self.assertEqual('Formalign.eu', brand_element.text)

        # She sees a form that says 'Paste in your alignment in FASTA format:'
        alignment_input = self.browser.find_element_by_css_selector('input#alignment')
        self.assertIsNotNone(self.browser.find_element_by_css_selector('label[for="alignment"]'))
        self.assertEqual(alignment_input.get_attribute('placeholder'), 'FASTA alignment')

        # She decides to give it a try. She pastes in her alignment of favorite proteins and submits it.
        with open(os.path.join(BASE_DIR, 'test_data/spa_align_clustal_omega.fasta'), 'r') as alignment_file:
            alignment_string = ''.join(line for line in alignment_file)
            alignment_input.send_keys(alignment_string)
        self.browser.find_element_by_css_selector('form button').click()
        self.fail('Incomplete Test')

        # She is redirected to a display page where she sees her alignment rendered in the default way.
