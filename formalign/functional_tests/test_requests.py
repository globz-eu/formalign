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

import requests
import time

__author__ = 'Stefan Dieterle'

start = time.time()
r = requests.get('http://localhost:8000/query-sequences/13/')
roundtrip = time.time() - start
print('13 query seq: ' + format(roundtrip))
start = time.time()
r = requests.get('http://localhost:8000/align-display/13/')
roundtrip = time.time() - start
print('13 alignment: ' + format(roundtrip))
start = time.time()
r = requests.get('http://localhost:8000/query-sequences/9/')
roundtrip = time.time() - start
print('9 query seq: ' + format(roundtrip))
start = time.time()
r = requests.get('http://localhost:8000/align-display/9/')
roundtrip = time.time() - start
print('9 alignment: ' + format(roundtrip))
