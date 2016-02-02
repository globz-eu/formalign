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
import time

from django.test import TestCase
from with_asserts.mixin import AssertHTMLMixin

from helper_funcs.helpers_bio import parse_fasta_alignment
from helper_funcs.helpers_test import file_to_string

from Bio.Alphabet.IUPAC import ExtendedIUPACProtein
from Bio.Alphabet import Gapped

from base.models import Alignment

__author__ = 'Stefan Dieterle'


class AlignDisplayTestCase(TestCase, AssertHTMLMixin):
    """
    Tests for alignment display
    """
    def setUp(self):
        """
        Creates a response from a GET request to /align-display/ with an alignment pk
        :param input_file: file containing alignment
        :return: response
        """
        name = 'A. tha. SPA family protein alignment'
        align_input = io.StringIO(file_to_string('ser_thr_kinase_family.fasta'))
        data = parse_fasta_alignment(align_input)
        for d in data:
            d.seq.alphabet = Gapped(ExtendedIUPACProtein())
        self.align = Alignment.objects.create_alignment(name, data)

    def test_align_display_speed(self):
        roundtrip = []
        for i in range(3):
            start = time.time()
            self.client.get('/align-display/' + str(self.align.id) + '/')
            elapsed = time.time() - start
            roundtrip.append(elapsed)
        average_time = sum(roundtrip) / len(roundtrip)
        self.assertTrue(average_time < 0.1, 'time: ' + format(average_time))
