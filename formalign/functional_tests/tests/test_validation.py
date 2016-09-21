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
import requests
from formalign.settings import SERVER_URL, TEST_CASE
from lxml import html
from io import StringIO
from base.forms import EMPTY_ERROR, FORMAT_ERROR, CHARACTER_ERROR, ALIGNMENT_ERROR, LESS_THAN_TWO_SEQS_ERROR

__author__ = 'Stefan Dieterle'


class InputValidationTestCase(TEST_CASE):
    """
    Tests input validation
    """
    def setUp(self):
        self.client = requests.Session()
        if SERVER_URL == 'liveserver':
            self.url = self.live_server_url
        else:
            self.url = SERVER_URL

    def tearDown(self):
        self.client.close()

    def alignment_validation(self, seqs):
        # She visits the Formalign.eu site
        r = self.client.get(self.url)
        csrftoken = self.client.cookies['csrftoken']
        index = html.parse(StringIO(r.text)).getroot()
        title = index.cssselect('title[id="head-title"]')
        brand = index.cssselect('a[class="navbar-brand"]')

        # page displays no error message
        self.assertEqual(r.status_code, 200, r.status_code)

        # User sees she's on the right page because she can see the name of the site in the title and the brand.
        self.assertEqual('Formalign.eu Home', title[0].text_content(), title[0].text_content())
        self.assertEqual(self.url + '/', r.url, r.url)
        self.assertEqual('Formalign.eu', brand[0].text_content(), brand[0].text_content())

        for t in ['protein', 'DNA']:
            for s in seqs:
                # she submits the invalid alignment
                alignment_string = file_to_string('%s%s.fasta' % (t, s['seq'])) if s['seq'] else None
                self.client.headers.update({'referer': self.url})
                r = self.client.post(self.url,
                                     data={'csrfmiddlewaretoken': csrftoken, 'seq_type': 'Protein',
                                           'align_input': alignment_string})
                index = html.parse(StringIO(r.text)).getroot()
                title = index.cssselect('title[id="head-title"]')
                brand = index.cssselect('a[class="navbar-brand"]')
                error_text = index.cssselect('ul[class="errorlist"]')[0].cssselect('li')[0].text_content()

                # response contains no error code
                self.assertEqual(r.status_code, 200, r.status_code)

                # user is redirected to the index page
                self.assertEqual('Formalign.eu Home', title[0].text_content(), title[0].text_content())
                self.assertEqual(self.url + '/', r.url, r.url)
                self.assertEqual('Formalign.eu', brand[0].text_content(), brand[0].text_content())

                # error text is displayed on index page
                self.assertEqual(s['error'], error_text, error_text)

    def test_alignment_validation(self):
        seqs = [
            {'seq': None, 'error': EMPTY_ERROR},
            {'seq': '_invalid_characters', 'error': '%ssequence1' % CHARACTER_ERROR},
            {'seq': '_too_few_sequences', 'error': LESS_THAN_TWO_SEQS_ERROR},
            {'seq': '_invalid_alignment', 'error': ALIGNMENT_ERROR}
        ]
        self.alignment_validation(seqs)

    def test_format_validation(self):
        seqs = [
            {'seq': '_invalid_fasta', 'error': FORMAT_ERROR},
        ]
        self.alignment_validation(seqs)
