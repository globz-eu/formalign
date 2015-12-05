from django.shortcuts import render
from base.forms import QueryForm
import os
from formalign.settings import BASE_DIR
from pprint import pprint

__author__ = 'Stefan Dieterle'


def index(request):
    """
    Serves home page
    :param request:
    :return:
    """
    form = QueryForm()
    return render(request, 'base/index.html', {'form': form})


def seq_display(request):
    """
    Serves query display page
    :param request:
    :return:
    """
    if request.method == 'POST':
        form = QueryForm(request.POST)
        if form.is_valid():
            return render(
                request,
                'base/query_display.html',
                {'query_seqs': parse_fasta(form.cleaned_data['align_input'].strip('\n'))}
            )
        else:
            return render(request, 'base/index.html', {'form': form})


def parse_fasta(fasta):
    """
    Parses fasta file into list of dicts of metadata and sequences
    :param fasta: fasta string
    :return: fasta_list [{'meta': '>sequence meta', 'seq': 'SEQUENCE'} ... ]
    """
    fasta_split = fasta.split('>')
    fasta_list = [f.split('\n') for f in fasta_split if f]
    # Remove some ghost whitespaces with: [''.join(fs.split()) for fs in f[1:]]
    fasta_dict_list = [
        {
           'meta': f[0].strip(),
           'seq':  ''.join([''.join(fs.split()) for fs in f[1:]])
        } for f in fasta_list
        ]
    return fasta_dict_list


if __name__ == '__main__':
    with open(os.path.join(BASE_DIR, 'test_data/short.fasta'), 'r') as alignment_file:
        input_seqs = alignment_file.read()
    parsed = parse_fasta('>Short sequence1\nMKERBGWAQ--Q\nGKKPWRF--EEW\n>Short sequence2\nMKERBGWA-SYQ\nGKKPWRFAQ-EW')
    print('parsed: ')
    pprint(parsed)
