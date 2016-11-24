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

__author__ = 'Stefan Dieterle'


def alignment_formatting(seqs_meta, tables):
    """
    tests that the alignments are displayed with the correct formatting
    :param seqs_meta: list of sequence meta data
    :param tables: lxml html object (align_table elements)
    """
    for nr, t in enumerate(tables):
        lines = t.find_class('al_ln')
        for i, l in enumerate(lines):
            assert 'seq_id' == l[0].attrib.get('class'), 'Unexpected class: %s' % l[0].attrib.get('class')
            assert seqs_meta[i] == l[0].text_content(), 'Unexpected meta data: %s' % l[0].text_content()
            assert len(l[1:]) <= 89, 'Line length was: %s' % len(l[1:])
            assert 'display_artifact' == l[-1].attrib.get('class'), \
                'Expected: display_artifact\nGot: %s' % l[-1].attrib.get('class')
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


def get_displayed_seqs(elems, alignment_length, cons=False, annot=False):
    """
    gets and recomposes sequences displayed on align-display page
    :return: re_seqs: recomposed sequences
    """
    # get displayed sequences
    seq_disp = []
    for els in elems:
        seq_disp_line = []
        for e in els.cssselect('td')[:-1]:
            if e.attrib.get('class') in ['residue S0', 'residue S1']:
                if annot:
                    seq_disp_line.append(e.attrib.get('class'))
                else:
                    seq_disp_line.append(e.text_content())
        if seq_disp_line:
            seq_disp.append(seq_disp_line)

    # recompose sequences
    if cons:
        cat_re_seq = []
        for j in range(alignment_length - 1, len(seq_disp), alignment_length):
            re_seq = [seq_disp[j] for j in range(alignment_length - 1, len(seq_disp), alignment_length)]
            cat_re_seq = []
            for r in re_seq:
                cat_re_seq.extend(r)
        return cat_re_seq
    else:
        re_seqs = []
        cat_re_seq = []
        if annot:
            length = alignment_length
        else:
            length = alignment_length + 1
        for i in range(0, length):
            for j in range(i, len(seq_disp), length):
                re_seq = [seq_disp[j] for j in range(i, len(seq_disp), length)]
                cat_re_seq = []
                for r in re_seq:
                    cat_re_seq.extend(r)
            re_seqs.append(cat_re_seq)
        return re_seqs
