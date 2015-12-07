import os
from formalign.settings import BASE_DIR

__author__ = 'Stefan Dieterle'


def file_to_string(file_name):
        with open(os.path.join(BASE_DIR, 'test_data/' + file_name), 'r') as input_seqs:
            alignment_string = input_seqs.read()
        return alignment_string
