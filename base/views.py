from django.shortcuts import render, redirect
from django.http import HttpResponseNotAllowed
from base.forms import QueryForm
from Bio.Align import AlignInfo
from base.models import save_alignment_to_db, get_multipleseqalignment_object_from_db

__author__ = 'Stefan Dieterle'


def index(request):
    """
    Serves home page
    :param request:
    :return:
    """
    if request.method == 'GET':
        form = QueryForm()
        return render(request, 'base/index.html', {'form': form})
    elif request.method == 'POST':
        form = QueryForm(request.POST)
        if form.is_valid():
            align = form.cleaned_data['align_input']
            pk = save_alignment_to_db('name', align)
            return redirect('/sequence-display/' + str(pk) + '/')
        else:
            return render(request, 'base/index.html', {'form': form})
    else:
        return HttpResponseNotAllowed(['POST', 'GET'])


def seq_display(request, align_id):
    """
    Serves query display page
    :param request:
    :param align_id:
    :return:
    """
    if request.method == 'GET':
        # split sequences in chunks of 80 characters
        length = 80
        align = get_multipleseqalignment_object_from_db(align_id)
        alphabets = {
            "Gapped(ExtendedIUPACProtein(), '-')": 'Protein',
            "Gapped(ExtendedIUPACDNA(), '-')": 'DNA',
        }
        alphabet = alphabets[str(align[0].seq.alphabet)]
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
            {'query_seqs': query_seqs, 'seq_type': alphabet, 'consensus': consensus}
        )

    else:
        return HttpResponseNotAllowed(['GET'])
