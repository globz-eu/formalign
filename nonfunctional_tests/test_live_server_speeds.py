from helper_funcs.helpers_test import file_to_string
from django.contrib.staticfiles.testing import StaticLiveServerTestCase
import requests
import time
from lxml import html
from io import StringIO
from statistics import mean


class BasicUserSpeedLiveServerTestCase(StaticLiveServerTestCase):
    """
    Tests basic user interaction with the app
    """

    def test_alignment_display_speeds(self):
        """
        Tests alignment display page
        """
        # User visits home page, submits a protein alignment and renders it
        files = [
            'spa_protein_alignment.fasta',
            'spa1_protein_alignment.fasta',
            'ser_thr_kinase_family.fasta'
        ]
        for f in files:
            q_display = []
            a_display = []
            for i in range(3):
                self.client = requests.Session()
                self.client.get(self.live_server_url)
                csrftoken = self.client.cookies['csrftoken']
                alignment_string = file_to_string(f)
                start = time.time()
                r = self.client.post(self.live_server_url,
                                     data={'csrfmiddlewaretoken': csrftoken, 'seq_type': 'protein',
                                           'cons_type': 'identity',
                                           'align_input': alignment_string, 'custom_data': 'custom'})
                roundtrip1 = time.time() - start
                q_display.append(roundtrip1)

                render_form = html.parse(StringIO(r.text)).getroot().cssselect('form[id="render"]')
                start = time.time()
                r = self.client.get(self.live_server_url + render_form[0].attrib.get('action'))
                roundtrip2 = time.time() - start
                a_display.append(roundtrip2)

                self.client.close()
            self.assertTrue(mean(q_display) <= 3, 'query sequences display of %s took: %s' % (f, mean(q_display)))
            self.assertTrue(mean(a_display) <= 3, 'alignment display of %s took: %s' % (f, mean(a_display)))
