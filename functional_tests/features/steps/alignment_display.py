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

from behave import then

from helper_funcs.helpers_test import file_to_string

__author__ = 'Stefan Dieterle'


@then(r'the alignment is displayed with 80 characters per line in blocks of 10 with sequence IDs')
def check_alignment_formatting(context):
    seqs_meta = file_to_string('spa_protein_alignment_meta.txt').splitlines()
    tables = context.display.find_class('align_table')
    for nr, t in enumerate(tables):
        lines = t.find_class('al_ln')
        for i, l in enumerate(lines):
            assert 'seq_id' == l[0].attrib.get('class'), 'Unexpected class: %s' % l[0].attrib.get('class')
            assert seqs_meta[i] == l[0].text_content(), 'Unexpected meta data: %s' % l[0].text_content()
            assert len(l[1:]) <= 89, 'Line length was: %s' % len(l[1:])
            assert 'display_artifact' == l[-1].attrib.get('class')
            if len(l[1:]) >= 12:
                for j in range(1, len(l[1:]) % 10):
                    assert 'block_sep' == l[j * 11].attrib.get('class'), \
                        'line number: %s, line length: %s, column: %s' % (i, len(l), format(j * 11))

                    for k in range(j * 11 - 10, j * 11):
                        assert l[k].attrib.get('class') == 'residue S0' or l[k].attrib.get('class') == 'residue S1', \
                            'class: %s, table: %s, line: %s, column: %s' % (l[k].attrib.get('class'), nr, i, k)
            else:
                for j in range(1, len(l)):
                    assert l[j].attrib.get('class') == 'residue S0' or l[j].attrib.get('class') == 'residue S1', \
                        l[j].attrib.get('class')
