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

import io
from django.test import TestCase
from django.core.urlresolvers import resolve
from helper_funcs.helpers_test import file_to_string
from helper_funcs.helpers_bio import parse_fasta_alignment
from Bio.Alphabet import Gapped
from Bio.Alphabet.IUPAC import ExtendedIUPACProtein
from base.views import index, seq_display, align_display
from base.models import Alignment

__author__ = 'Stefan Dieterle'


class BaseURLsTestCase(TestCase):
    """
    Tests URLs resolve to the right views
    """
    def setUp(self):
        name = 'A. tha. SPA family protein alignment'
        align_input = io.StringIO(file_to_string('spa_protein_alignment.fasta'))
        data = parse_fasta_alignment(align_input)
        for d in data:
            d.seq.alphabet = Gapped(ExtendedIUPACProtein())
        align = Alignment.objects.create_alignment(name, data)
        self.slug = align.slug

    def test_root_url_uses_index_view(self):
        """
        Tests that the root of the site resolves to the correct view function
        """
        root = resolve('/')
        self.assertEqual(root.func, index)

    def test_query_seq_url_uses_seq_display_view(self):
        """
        Test that the /query-sequences/ URL resolves to the correct view function
        """
        query_seq = resolve('/query-sequences/%s' % self.slug)
        self.assertEqual(query_seq.func, seq_display)

    def test_align_display_url_uses_align_display_view(self):
        """
        Test that the /align-display/ URL resolves to the correct view function
        """
        align_disp = resolve('/align-display/%s' % self.slug)
        self.assertEqual(align_disp.func, align_display)
