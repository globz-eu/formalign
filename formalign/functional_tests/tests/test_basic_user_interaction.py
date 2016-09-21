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

__author__ = 'Stefan Dieterle'


class BasicUserTestCase(TEST_CASE):
    """
    Tests basic user interaction with the app
    """

    def setUp(self):
        self.client = requests.Session()
        if SERVER_URL == 'liveserver':
            self.url = self.live_server_url
        else:
            self.url = SERVER_URL

    def tearDown(self):
        self.client.close()

    def test_index(self):
        """
        Tests index page
        """
        # Lambda user is a biologist who has to make a nice figure containing a multiple alignment for a presentation.
        # She visits the formalign.eu site.
        r = self.client.get(self.url)
        csrftoken = self.client.cookies['csrftoken']
        index = html.parse(StringIO(r.text)).getroot()
        title = index.cssselect('title[id="head-title"]')
        brand = index.cssselect('a[class="navbar-brand"]')
        align_input = index.cssselect('textarea[id="id_align_input"]')
        prot_button = index.cssselect('input[id="id_seq_type_0"]')
        prot_button_label = index.cssselect('label[for="id_seq_type_0"]')
        dna_button = index.cssselect('input[id="id_seq_type_1"]')

        # page displays no error message
        self.assertEqual(r.status_code, 200, r.status_code)

        # User sees she's on the right page because she can see the name of the site in the title and the brand.
        self.assertEqual('Formalign.eu Home', title[0].text_content(), title[0].text_content())
        self.assertEqual(self.url + '/', r.url, r.url)
        self.assertEqual('Formalign.eu', brand[0].text_content(), brand[0].text_content())

        # She sees a form that says 'Paste in your alignment'
        self.assertEqual('Paste in your alignment:(FASTA, clustalw, stockholm or phylip)',
                         align_input[0].label.text_content(), align_input[0].label.text_content())

        # There's a textarea with a placeholder saying 'Alignment'
        self.assertEqual('Alignment (FASTA, clustalw, stockholm or phylip)',
                         align_input[0].attrib.get('placeholder'), align_input[0].attrib.get('placeholder'))

        # She sees two radio buttons for DNA and protein
        self.assertEqual('Input sequence type:', prot_button[0].label.text_content(),
                         prot_button[0].label.text_content())
        self.assertEqual(' Protein', prot_button_label[1].text_content(), prot_button_label[1].text_content())
        self.assertEqual(' DNA', dna_button[0].label.text_content(), dna_button[0].label.text_content())

        # She sees that the DNA button is selected by default
        self.assertTrue(prot_button[0].checkable, 'Protein button is not checkable')
        self.assertFalse(prot_button[0].checked, 'Protein button was checked by default')
        self.assertTrue(dna_button[0].checkable, 'DNA button is not checkable')
        self.assertTrue(dna_button[0].checked, 'DNA button was not checked by default')

        # She checks the protein button and pastes in an alignment
        alignment_string = file_to_string('spa_protein_alignment.fasta')
        r = self.client.post(self.url,
                             data={'csrfmiddlewaretoken': csrftoken, 'seq_type': 'Protein',
                                   'align_input': alignment_string})
        self.assertEqual(r.status_code, 200, r.status_code)
        display = html.parse(StringIO(r.text)).getroot()
        title = display.cssselect('title[id="head-title"]')

        # She is redirected to the sequence display page
        self.assertEqual('Formalign.eu Sequence Display', title[0].text_content(), title[0].text_content())
        self.assertEqual('%s/query-sequences' % self.url, '/'.join(r.url.split('/')[:-2]),
                         '/'.join(r.url.split('/')[:-2]))

    def test_query_seqs_display(self):
        """
        Tests the query sequences display page
        """
        # User visits home page and submits a protein alignment
        self.client.get(self.url)
        csrftoken = self.client.cookies['csrftoken']
        alignment_string = file_to_string('spa_protein_alignment.fasta')
        r = self.client.post(self.url,
                             data={'csrfmiddlewaretoken': csrftoken, 'seq_type': 'Protein',
                                   'align_input': alignment_string})
        display = html.parse(StringIO(r.text)).getroot()
        title = display.cssselect('title[id="head-title"]')
        sequence_lines = display.cssselect('p[class="query_seq_display"]')
        sequence_meta = display.find_class('query_seq_meta')
        render_form = display.cssselect('form[id="render"]')

        # She is redirected to a page showing the submitted sequences from her alignment and a simple consensus sequence
        self.assertEqual(r.status_code, 200, r.status_code)
        self.assertEqual('%s/query-sequences' % self.url, '/'.join(r.url.split('/')[:-2]),
                         '/'.join(r.url.split('/')[:-2]))
        self.assertEqual('Formalign.eu Sequence Display', title[0].text_content(), title[0].text_content())
        self.assertIsNotNone(sequence_lines, 'sequences are empty')
        for l in sequence_lines:
            self.assertTrue(len(l.text_content()) <= 80, l.text_content())
        seqs = file_to_string('spa_protein_alignment_display.txt').splitlines()
        for i, a in enumerate(seqs):
            self.assertEqual(a, sequence_lines[i].text_content(),
                             '%s: %s' % (format(i), sequence_lines[i].text_content()))
        seqs_meta = file_to_string('spa_protein_alignment_display_meta.txt').splitlines()
        for i, a in enumerate(seqs_meta):
            self.assertEqual(a, sequence_meta[i].text_content(), sequence_meta[i].text_content())

        # She is happy with the result, sees a "Render" button and clicks it.
        self.assertEqual('get', render_form[0].attrib.get('method'), render_form[0].attrib.get('method'))
        self.assertEqual('align-display', render_form[0].attrib.get('action').split('/')[1],
                         render_form[0].attrib.get('action').split('/')[1])
        self.assertTrue(render_form[0].attrib.get('action').split('/')[-2].isdigit(),
                        render_form[0].attrib.get('action').split('/')[-2])
        r = self.client.get(self.url + render_form[0].attrib.get('action'))

        # She is redirected to the render page
        self.assertEqual(r.status_code, 200, r.status_code)
        align = html.parse(StringIO(r.text)).getroot()
        title = align.cssselect('title[id="head-title"]')
        self.assertEqual('Formalign.eu Alignment Display', title[0].text_content(), title[0].text_content())
        self.assertEqual('%s/align-display' % self.url, '/'.join(r.url.split('/')[:-2]),
                         '/'.join(r.url.split('/')[:-2]))

    def test_alignment_display(self):
        """
        Tests alignment display page
        """
        # User visits home page, submits a protein alignment and renders it
        self.client.get(self.url)
        csrftoken = self.client.cookies['csrftoken']
        alignment_string = file_to_string('spa_protein_alignment.fasta')
        r = self.client.post(self.url,
                             data={'csrfmiddlewaretoken': csrftoken, 'seq_type': 'Protein',
                                   'align_input': alignment_string})
        render_form = html.parse(StringIO(r.text)).getroot().cssselect('form[id="render"]')
        r = self.client.get(self.url + render_form[0].attrib.get('action'))

        # She is redirected to the render page
        self.assertEqual(r.status_code, 200, r.status_code)
        align = html.parse(StringIO(r.text)).getroot()
        title = align.cssselect('title[id="head-title"]')
        self.assertEqual('Formalign.eu Alignment Display', title[0].text_content(), title[0].text_content())
        self.assertEqual('%s/align-display' % self.url, '/'.join(r.url.split('/')[:-2]),
                         '/'.join(r.url.split('/')[:-2]))

        # She sees the alignment displayed with 80 characters per line in blocks of 10 with sequence ids
        seqs_meta = file_to_string('spa_protein_alignment_meta.txt').splitlines()
        tables = align.find_class('align_table')
        for nr, t in enumerate(tables):
            lines = t.find_class('al_ln')
            for i, l in enumerate(lines):
                self.assertEqual('seq_id', l[0].attrib.get('class'), l[0].attrib.get('class'))
                self.assertEqual(seqs_meta[i], l[0].text_content(), l[0].text_content())
                self.assertTrue(len(l[1:]) <= 89, len(l[1:]))
                self.assertEqual('display_artifact', l[-1].attrib.get('class'))
                if len(l[1:]) >= 12:
                    for j in range(1, len(l[1:]) % 10):
                        self.assertEqual('block_sep', l[j * 11].attrib.get('class'),
                                         'line number: %s, line length: %s, column: %s' % (i, len(l), format(j * 11)))

                        for k in range(j * 11 - 10, j * 11):
                            self.assertTrue(
                                l[k].attrib.get('class') == 'residue S0' or l[k].attrib.get('class') == 'residue S1',
                                'class: %s, table: %s, line: %s, column: %s' % (l[k].attrib.get('class'), nr, i, k))
                else:
                    for j in range(1, len(l)):
                        self.assertTrue(
                            l[j].attrib.get('class') == 'residue S0' or l[j].attrib.get('class') == 'residue S1',
                            l[j].attrib.get('class'))

        # She is quite happy with the result and decides to navigate back to the home page
        home_button = align.cssselect('a[class="navbar-brand"]')
        self.assertEqual('Formalign.eu', home_button[0].text_content(), home_button[0].text_content())
        r = self.client.get(self.url + home_button[0].attrib.get('href'))
        self.assertEqual(r.status_code, 200, r.status_code)
        home = html.parse(StringIO(r.text)).getroot()
        title = home.cssselect('title[id="head-title"]')

        # User sees she's that she's back on the home page she can see the name of the site in the title and the brand.
        self.assertEqual('Formalign.eu Home', title[0].text_content(), title[0].text_content())
        self.assertEqual(self.url + '/', r.url, r.url)
