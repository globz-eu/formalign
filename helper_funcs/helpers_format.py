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
    :return: formatted line split object:
        for split_type = sequence:  [
                                        {'seq': [
                                                    SeqRecord(sequence_1, 1st to line_length residues),
                                                    SeqRecord(sequence_1, line_length+1 to 2*line_length residues),
                                                    ...,
                                                    SeqRecord(sequence_1, left over residues),
                                                ]
                                        },
                                        ...,
                                        {'seq': [
                                                    SeqRecord(sequence_n, 1st to line_length residues),
                                                    SeqRecord(sequence_n, line_length+1 to 2*line_length residues),
                                                    ...,
                                                    SeqRecord(sequence_n, left over residues),
                                                ]
                                        },
                                        <if present in alignment: {'seq': [
                                                    SeqRecord(consensus, 1st to line_length residues),
                                                    SeqRecord(consensus, line_length+1 to 2*line_length residues),
                                                    ...,
                                                    SeqRecord(consensus, left over residues),
                                                ]
                                        }>
                                    ]

        for split_type = alignment: [
                                        [
                                            SeqRecord(sequence_1, 1st to line_length residues),
                                            ...,
                                            SeqRecord(sequence_n, 1st to line_length residues),
                                            <if present in alignment: SeqRecord(consensus, 1st to line_length residues)>
                                        ],
                                        [
                                            SeqRecord(sequence_1, line_length+1 to 2*line_length residues),
                                            ...,
                                            SeqRecord(sequence_n, line_length+1 to 2*line_length residues),
                                            <if present in alignment:
                                            SeqRecord(consensus, line_length+1 to 2*line_length residues)>
                                        ],
                                        ...,
                                        [
                                            SeqRecord(sequence_1, left over residues),
                                            ...,
                                            SeqRecord(sequence_n, left over residues),
                                            <if present in alignment: SeqRecord(consensus, left over residues)>
                                        ],
                                    ]
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


def split_lines_in_blocks(alignment, block_size=10):
    """
    splits lines of an alignment in blocks of block_size residues
    :param: alignment: formatted line split object returned by split_lines
    :return: formatted block split object:
    [
        [                                                                       # first sequence block
            [                                                                   # first line of first sequence
                SeqRecord.id(sequence1),
                [
                    [                                                           # first block_size residues of
                                                                                # sequence 1 as zip object
                        (
                            res(1st) annotation,
                            res(1st)
                        ),
                        ...,
                        (
                            res(block_size th),
                            res(block_size th) annotation
                        )
                    ],
                    ...,
                    [                                                           # last block_size residues
                                                                                # as zip object
                        (
                            res(line length th - block_size th) annotation,
                            res(line length th - block_size th)
                        ),
                        ...,
                        (
                            res(line length th),
                            res(line length th) annotation
                        )
                    ],
                ],
            ],
            ...,
            [                                                                   # first line of last sequence
                SeqRecord.id(sequence_n),                                       # or of consensus sequence
                    [                                                           # first block_size residues of
                                                                                # sequence n as zip object
                        (
                            res(1st) annotation,
                            res(1st)
                        ),
                        ...,
                        (
                            res(block_size th),
                            res(block_size th) annotation
                        )
                    ],
                    ...,
                    [                                                           # last block_size residues
                                                                                # as zip object
                        (
                            res(line length th - block_size th) annotation,
                            res(line length th - block_size th)
                        ),
                        ...,
                        (
                            res(line length th),
                            res(line length th) annotation
                        )
                    ],
                ],
            ],
        ],
        ...,
        [
            [                                                                   # last line of first sequence
                SeqRecord.id(sequence1),
                [
                    [],                                                         # first block_size residues of sequence 1
                    ...,                                                        # or remaining residues
                    []
                ],
            ],
            ...,
            [                                                                   # last line of last sequence
                SeqRecord.id(sequence1),                                        # or of consensus sequence
                [
                    [],                                                         # first block_size residues of sequence 1
                    ...,                                                        # or remaining residues
                    []
                ]
            ]
        ]
    ]
    """
    seqs_blocks = [
        [
            [
                ls.id, [
                    zip(
                        ls.letter_annotations['eq'][j:j + block_size],
                        list(ls[j:j + block_size])
                    ) for j in range(0, len(ls), block_size)
                ]
            ] for ls in ln] for ln in alignment
        ]
    return seqs_blocks


def annotate(alignment):
    """
    Annotates alignment according to consensus sequence. Black if equal to consensus or white if not equal for
    sequences. White if X or - (undefined residues), else Black (defined residues)
    :param alignment: MultipleSeqAlignment object with consensus
    :return: MultipleSeqAlignment object with consensus and letter_annotations
    """
    cons_seq_annot = {'eq': [0 if c in ['-', 'X'] else 1 for c in alignment[-1].seq]}
    alignment[-1].letter_annotations = cons_seq_annot
    for a in alignment[:-1]:
        a.letter_annotations = {
            'eq': [1 if (a[i] == alignment[-1][i] and alignment[-1][i] not in ['-', 'X'])
                   else 0
                   for i in range(0, len(a))]
        }
    return alignment
