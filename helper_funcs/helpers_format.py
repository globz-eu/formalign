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
