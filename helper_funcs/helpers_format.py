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


def split_lines(alignment, line_length=80, split_type='sequence'):
    """
    Splits sequences in lines of line_length characters or alignment in alignment blocks of line_length characters
    :param alignment: MultipleSeqAlignment object
    :param line_length: length of lines
    :param split_type: 'sequence' for sequence split, 'alignment' for alignment split
    :return: formatted split object
    """
    seq_lines = []
    if split_type == 'sequence':
        seq_lines = [
                        {
                            'description': f.description,
                            'id': f.id,
                            'name': f.name,
                            'seq': [
                                        f.seq[i:i + line_length] for i in range(0, len(f.seq), line_length)
                                    ]
                        } for f in alignment
                    ]

    elif split_type == 'alignment':
        seq_lines = [
                        [
                            f[i:i + line_length] for f in alignment
                        ] for i in range(0, len(alignment[0].seq), line_length)
                    ]
    return seq_lines


def split_lines_in_blocks():
    pass


def consensus_annotate():
    pass
