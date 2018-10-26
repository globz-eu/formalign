import os
import io
from formalign.settings import BASE_DIR
from helper_funcs.bio.helpers import parse_fasta_alignment, consensus_add
from helper_funcs.helpers_format import annotate


def file_to_string(file_name):
    """
    converts test file contents to a string
    :param file_name: file to convert to a string
    :return: string
    """
    with open(os.path.join(BASE_DIR, 'test_data/' + file_name), 'r') as input_seqs:
        alignment_string = input_seqs.read()
    return alignment_string


def sequences_file(alignment_file, file_name='default'):
    """
    creates a file with one sequence per line from a fasta file containing an
    alignment
    :param file_name: output file name
    :param alignment_file: alignments file name
    :return: 0 if successful
    """
    if file_name == 'default':
        file_name = alignment_file.replace('.fasta', '_seqs.txt')
    alignment_string = file_to_string(alignment_file)
    alignment = parse_fasta_alignment(io.StringIO(alignment_string))
    alignment_list = [''.join(list(al.seq)) + '\n' for al in alignment]
    with open(os.path.join(BASE_DIR, 'test_data/' + file_name), 'w') as align:
        for a in alignment_list:
            align.write(a)
    return 0


def sequences_annotations_file(alignment_file, file_name='default'):
    """
    creates a file with an annotation sequence of a sequence per line from a
    fasta file containing an alignment
    :param file_name: output file name
    :param alignment_file: alignment file name
    :return: 0 if successful
    """
    if file_name == 'default':
        file_name = alignment_file.replace('.fasta', '_seqs_annot.txt')
    alignment_string = file_to_string(alignment_file)
    alignment = parse_fasta_alignment(io.StringIO(alignment_string))
    alignment = consensus_add(alignment)
    alignment = annotate(alignment)
    alignment_list = [''.join(map(str, al.letter_annotations['eq'])) + '\n' for al in alignment]
    with open(os.path.join(BASE_DIR, 'test_data/' + file_name), 'w') as align:
        for a in alignment_list:
            align.write(a)
    return 0


def consensus_file(alignment_file, file_name='default'):
    """
    creates a file with the consensus sequence for the alignment
    :param file_name: output file name
    :param alignment_file: alignment file name
    :return: 0 if successful
    """
    if file_name == 'default':
        file_name = alignment_file.replace('.fasta', '_consens.txt')
    alignment_string = file_to_string(alignment_file)
    alignment = parse_fasta_alignment(io.StringIO(alignment_string))
    alignment = consensus_add(alignment)
    alignment_list = [''.join(list(al.seq)) + '\n' for al in alignment]
    with open(os.path.join(BASE_DIR, 'test_data/' + file_name), 'w') as align:
        for a in alignment_list:
            align.write(a)
    return 0
