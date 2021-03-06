import time
from formalign.settings import CHROME_DRIVER, SERVER_URL, BROWSER_DELAY, TEST_CASE, FIREFOX_BINARY
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.desired_capabilities import DesiredCapabilities
import pyperclip
from helper_funcs.helpers_test import file_to_string


class BasicUserTestCaseChrome(TEST_CASE):
    def setUp(self):
        self.browser = webdriver.Chrome(
            CHROME_DRIVER
        )
        self.delay = float(BROWSER_DELAY)
        if SERVER_URL == 'liveserver':
            self.url = self.live_server_url
            self.sleep = self.delay * 5
        else:
            self.url = SERVER_URL
            self.sleep = self.delay

    def tearDown(self):
        self.browser.quit()

    def test_demo_alignment(self):
        """
        Tests submission of the demo alignment
        """
        # Lambda user is a biologist who has to make a nice figure containing a multiple alignment for a presentation.
        # She visits the formalign.eu site.
        self.browser.get(self.url + '/')

        # User sees she's on the right page because she can see the name of the site in the heading.
        self.assertEqual(self.browser.title, 'Formalign.eu Home', self.browser.title)

        # She sees a use demo alignment button and clicks it
        submit_demo = self.browser.find_element_by_id('submit-demo')
        self.assertEqual('custom_data', submit_demo.get_attribute('name'), submit_demo.get_attribute('name'))
        self.assertEqual('demo', submit_demo.get_attribute('value'), submit_demo.get_attribute('value'))
        submit_demo.click()

        # Wait for browser
        time.sleep(self.sleep * 10)

        # She is redirected to a page showing the submitted sequences from her alignment and a simple consensus sequence
        self.assertEqual(self.browser.title, 'Formalign.eu Sequence Display', self.browser.title)
        seq_content = self.browser.find_elements_by_css_selector('.query_seq_display')
        self.assertIsNotNone(seq_content)
        # for f in seq_content:
        #     self.assertTrue(len(f.text) <= 80)

        first_seq_info = self.browser.find_elements_by_css_selector('.query_seq_meta')[0]
        self.assertEqual(
            'DMD401_1-640:',
            first_seq_info.text,
            first_seq_info.text
        )
        first_seq_content = self.browser.find_elements_by_css_selector('.query_seq_display')[0]
        self.assertIsNotNone(first_seq_content)
        self.assertEqual(
            'LQLDTVLGEGEFGQVLKGFATEIAG---------LPGITTVAVKMLKKGSNSV------------EYMALLSEFQLLQEV',
            first_seq_content.text
        )

        consensus_seq = self.browser.find_elements_by_xpath(
            '//div[@class="query_seq bg-color-body"]'
        )[-1].find_elements_by_xpath('./p[@class="query_seq_display"]')[0]
        self.assertIsNotNone(consensus_seq)
        cons_seq = file_to_string('consensus_ser_thr_kin.txt')
        self.assertEqual(consensus_seq.text, cons_seq[:80])
        consensus_meta = self.browser.find_elements_by_xpath('//h3[@class="query_seq_meta bg-color-body"]')[-1]
        self.assertEqual(consensus_meta.text, 'consensus 70%:')

        # She is happy with the result, sees a "Render" button and clicks it.
        render_button = self.browser.find_element_by_css_selector('button#render-align')
        self.assertIsNotNone(render_button)
        render_button.click()
        # Wait for Firefox
        time.sleep(self.sleep * 10)

        # She is redirected to the alignment display page
        self.assertEqual('Formalign.eu Alignment Display', self.browser.title, self.browser.title)

        # She is quite happy with the result and decides to try with another alignment so she navigates back to the
        # home page
        home_button = self.browser.find_element_by_css_selector('.navbar-brand')
        home_button.click()
        # Wait for Firefox
        time.sleep(self.sleep * 2)

        self.assertEqual('Formalign.eu Home', self.browser.title, self.browser.title)

    def test_basic_user_experience(self):
        """
        Tests basic user interaction with formalign.eu site
        """
        # Lambda user is a biologist who has to make a nice figure containing a multiple alignment for a presentation.
        # She visits the formalign.eu site.
        self.browser.get(self.url + '/')

        # User sees she's on the right page because she can see the name of the site in the heading.
        self.assertEqual(self.browser.title, 'Formalign.eu Home', self.browser.title)
        brand_element = self.browser.find_element_by_css_selector('.navbar-brand')
        self.assertEqual('Formalign.eu', brand_element.text)

        # She sees a form that says 'Paste in your alignment in FASTA format:'
        alignment_input = self.browser.find_element_by_css_selector('textarea#id_align_input')
        self.assertIsNotNone(self.browser.find_element_by_css_selector('label[for="id_align_input"]'))
        self.assertEqual(
            'Alignment (FASTA, clustalw, stockholm or phylip)',
            alignment_input.get_attribute('placeholder'),
        )

        # She sees two radio buttons for DNA and protein
        dna_button = self.browser.find_element_by_css_selector('input#id_seq_type_1')
        self.assertIsNotNone(dna_button)
        protein_button = self.browser.find_element_by_css_selector('input#id_seq_type_0')
        self.assertIsNotNone(protein_button)

        # She sees that the DNA button is selected by default
        self.assertEqual(dna_button.is_selected(), True)

        # She clicks the Protein radio button and sees that it gets selected and the DNA button gets unselected
        protein_button.click()
        # Wait for Firefox
        time.sleep(self.sleep)

        self.assertEqual(protein_button.is_selected(), True)
        self.assertEqual(dna_button.is_selected(), False)

        # She sees two radio buttons for identity and substitution matrix
        # consensus
        identity_button = self.browser.find_element_by_css_selector('input#id_cons_type_0')
        self.assertIsNotNone(identity_button)
        subs_matrix_button = self.browser.find_element_by_css_selector('input#id_cons_type_1')
        self.assertIsNotNone(subs_matrix_button)

        # She sees that the identity button is selected by default
        self.assertEqual(identity_button.is_selected(), True)

        # She clicks the substitution matrix radio button and sees that it gets
        # selected and the identity button gets unselected
        subs_matrix_button.click()
        # Wait for Firefox
        time.sleep(self.sleep)

        self.assertEqual(subs_matrix_button.is_selected(), True)
        self.assertEqual(identity_button.is_selected(), False)

        # She goes back to identity
        identity_button.click()
        # Wait for Firefox
        time.sleep(self.sleep)

        self.assertEqual(subs_matrix_button.is_selected(), False)
        self.assertEqual(identity_button.is_selected(), True)

        # She pastes in a protein alignment to see what happens
        alignment_string = file_to_string('spa_protein_alignment.fasta')
        pyperclip.copy(alignment_string)
        alignment_input = self.browser.find_element_by_css_selector('textarea#id_align_input')
        alignment_input.send_keys(Keys.CONTROL, 'v')
        self.browser.find_element_by_id('submit-align').click()

        # Wait for Firefox
        time.sleep(self.sleep * 10)

        # She is redirected to a page showing the submitted sequences from her alignment and a simple consensus sequence
        self.assertEqual(self.browser.title, 'Formalign.eu Sequence Display', self.browser.title)
        seq_content = self.browser.find_elements_by_css_selector('.query_seq_display')
        self.assertIsNotNone(seq_content)
        for f in seq_content:
            self.assertTrue(len(f.text) <= 80)

        first_seq_info = self.browser.find_elements_by_css_selector('.query_seq_meta')[0]
        self.assertEqual(
            'NP_175717 NP_175717.1 SPA1-related 4 protein [Arabidopsis thaliana].:',
            first_seq_info.text,
            first_seq_info.text
        )
        first_seq_content = self.browser.find_elements_by_css_selector('.query_seq_display')[0]
        self.assertIsNotNone(first_seq_content)
        self.assertEqual(first_seq_content.text, '-' * 80)

        consensus_seq = self.browser.find_elements_by_xpath(
            '//div[@class="query_seq bg-color-body"]'
        )[-1].find_elements_by_xpath('./p[@class="query_seq_display"]')[0]
        self.assertIsNotNone(consensus_seq)
        cons_seq = file_to_string('consensus.txt')
        self.assertEqual(consensus_seq.text, cons_seq[:80])
        consensus_meta = self.browser.find_elements_by_xpath('//h3[@class="query_seq_meta bg-color-body"]')[-1]
        self.assertEqual(consensus_meta.text, 'consensus 70%:')

        # She is happy with the result, sees a "Render" button and clicks it.
        render_button = self.browser.find_element_by_css_selector('button#render-align')
        self.assertIsNotNone(render_button)
        render_button.click()
        # Wait for Firefox
        time.sleep(self.sleep * 10)

        # She is redirected to the alignment display page
        self.assertEqual('Formalign.eu Alignment Display', self.browser.title, self.browser.title)

        # She sees the alignment displayed with 80 characters per line in blocks of 10 with sequence ids
        s0 = self.browser.find_elements_by_xpath(
            '//tr[@class="al_ln"]'
        )[10].find_elements_by_xpath('./td[@class="residue S0"]')
        s1 = self.browser.find_elements_by_xpath(
            '//tr[@class="al_ln"]'
        )[10].find_elements_by_xpath('./td[@class="residue S1"]')
        self.assertEqual(len(s0) + len(s1), 80)
        sep = self.browser.find_elements_by_xpath(
            '//tr[@class="al_ln"]'
        )[10].find_elements_by_xpath('./td[@class="block_sep"]')
        self.assertEqual(len(sep), 8)

        # She is quite happy with the result and decides to try with another alignment so she navigates back to the
        # home page
        home_button = self.browser.find_element_by_css_selector('.navbar-brand')
        home_button.click()
        # Wait for Firefox
        time.sleep(self.sleep * 2)

        self.assertEqual('Formalign.eu Home', self.browser.title, self.browser.title)

        # She wants to upload a protein stockholm alignment this time from a file
        # She clicks the Protein radio button and sees that it gets selected and the DNA button gets unselected
        protein_button = self.browser.find_element_by_css_selector('input#id_seq_type_0')
        protein_button.click()
        self.assertEqual(protein_button.is_selected(), True)


class BasicUserTestCaseFirefox(BasicUserTestCaseChrome):
    def setUp(self):
        caps = DesiredCapabilities.FIREFOX
        caps['marionette'] = True
        caps['binary'] = FIREFOX_BINARY
        self.browser = webdriver.Firefox(capabilities=caps)
        self.delay = float(BROWSER_DELAY)
        if SERVER_URL == 'liveserver':
            self.url = self.live_server_url
            self.sleep = self.delay * 10
        else:
            self.url = SERVER_URL
            self.sleep = self.delay * 2.5

    def tearDown(self):
        self.browser.quit()
