"""
=====================================================================
Django app deployment scripts
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

from django.contrib.staticfiles.testing import StaticLiveServerTestCase
from configuration import CHROME_DRIVER
from selenium import webdriver

__author__ = 'Stefan Dieterle'


class BasicUserTestCase(StaticLiveServerTestCase):
    """
    Tests basic user interaction with the app
    """

    def setUp(self):
        self.browser = webdriver.Chrome(
            CHROME_DRIVER
        )
        # self.browser = webdriver.Firefox()
        self.browser.implicitly_wait(2)

    def tearDown(self):
        # time.sleep(5)
        self.browser.quit()

    def test_file_upload(self):
        # User visits the formalign.eu site
        self.browser.get(self.live_server_url + '/')

        # She wants to upload a protein stockholm alignment this time from a file
        # She clicks the Protein radio button and sees that it gets selected and the DNA button gets unselected
        protein_button = self.browser.find_element_by_css_selector('input#id_seq_type_0')
        protein_button.click()
        self.assertEqual(protein_button.is_selected(), True)

        # She sees a file upload button
        alignment_input = self.browser.find_element_by_css_selector('file_upload#id_align_input')
        alignment_input.click()

        # She browses to her file and selects it

        # She submits her file
        self.browser.find_element_by_id('submit-align').click()
        self.fail('Incomplete Test')
