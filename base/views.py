from django.shortcuts import render
from base.forms import QueryForm
from Bio.Align import AlignInfo

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
            # split sequences in chunks of 80 characters
            length = 80
            align = form.cleaned_data['align_input']
            query_seqs = [
                {
                    'meta': f.description,
                    'seq': [
                        f.seq[i:i+length] for i in range(0, len(f.seq), length)
                        ]
                } for f in align
                ]
            cons_seq = AlignInfo.SummaryInfo(align).gap_consensus()
            consensus = {'meta': 'Consensus', 'seq': [
                        cons_seq[i:i+length] for i in range(0, len(cons_seq), length)
                        ]}
            return render(
                request,
                'base/query_display.html',
                {'query_seqs': query_seqs, 'seq_type': form.cleaned_data['seq_type'], 'consensus': consensus}
            )
        else:
            return render(request, 'base/index.html', {'form': form})
