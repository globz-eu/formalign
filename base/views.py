from django.shortcuts import render
from base.forms import QueryForm

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
                {'query_seqs': parse_fasta(form.cleaned_data['align_input'])}
            )
        else:
            return render(request, 'base/index.html', {'form': form})


def parse_fasta(fasta):
    """
    Parses fasta file into list of dicts of metadata and sequences
    :param fasta: fasta string
    :return: fasta_list [{'meta': '>sequence meta', 'seq': 'SEQUENCE'} ... ]
    """
    fasta_split = fasta.split('\n')
    fasta_list = [{'meta': fasta_split[n], 'seq': fasta_split[n+1]} for n in range(0, len(fasta_split), 2)]
    return fasta_list
