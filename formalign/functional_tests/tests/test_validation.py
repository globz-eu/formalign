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

from django.contrib.staticfiles.testing import StaticLiveServerTestCase
from acceptance_tests.requests_tests.test_validation import InputValidationTestCase
import requests

__author__ = 'Stefan Dieterle'


class InputValidationTestCaseLiveServer(InputValidationTestCase, StaticLiveServerTestCase):
    """
    Tests input validation
    """
    def setUp(self):
        self.client = requests.Session()
        self.url = self.live_server_url

    def tearDown(self):
        self.client.close()
