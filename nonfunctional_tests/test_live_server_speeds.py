"""
=====================================================================
Formalign.eu format and display multiple sequence alignments
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

from helper_funcs.helpers_test import file_to_string
from django.contrib.staticfiles.testing import StaticLiveServerTestCase
import requests
import time
from lxml import html
from io import StringIO
from statistics import mean

__author__ = 'Stefan Dieterle'


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
                                     data={'csrfmiddlewaretoken': csrftoken, 'seq_type': 'Protein',
                                           'align_input': alignment_string, 'custom_data': 'custom'})
                roundtrip1 = time.time() - start
                q_display.append(roundtrip1)

                render_form = html.parse(StringIO(r.text)).getroot().cssselect('form[id="render"]')
                start = time.time()
                r = self.client.get(self.live_server_url + render_form[0].attrib.get('action'))
                roundtrip2 = time.time() - start
                a_display.append(roundtrip2)

                self.client.close()
            self.assertTrue(mean(q_display) <= 2, 'query sequences display of %s took: %s' % (f, mean(q_display)))
            self.assertTrue(mean(a_display) <= 2, 'alignment display of %s took: %s' % (f, mean(a_display)))
