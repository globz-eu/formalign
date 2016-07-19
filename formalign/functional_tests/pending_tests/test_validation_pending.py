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

from formalign.functional_tests.tests.test_validation import InputValidationTestCase
from base.forms import FORMAT_ERROR

__author__ = 'Stefan Dieterle'


class InputValidationTestCasePending(InputValidationTestCase):
    def test_format_validation(self):
        seqs = [
            {'seq': '_invalid_fasta', 'error': FORMAT_ERROR},
            {'seq': '_invalid_clustal', 'error': FORMAT_ERROR},
            {'seq': '_invalid_ig', 'error': FORMAT_ERROR},
            {'seq': '_invalid_nexus', 'error': FORMAT_ERROR},
            {'seq': '_invalid_phylip', 'error': FORMAT_ERROR},
            {'seq': '_invalid_stockholm', 'error': FORMAT_ERROR},
        ]
        self.alignment_validation(seqs)
