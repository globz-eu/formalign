from django.shortcuts import render, redirect
from django.http import HttpResponseNotAllowed
from base.forms import QueryForm
from Bio.Align import AlignInfo
from Bio.SeqRecord import SeqRecord
from base.models import Alignment

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
            save_align = Alignment.objects.create_alignment('name', align)
            pk = save_align.id
            return redirect('/query-sequences/' + str(pk) + '/')
        else:
            return render(request, 'base/index.html', {'form': form})
    else:
        return HttpResponseNotAllowed(['POST', 'GET'])


def seq_display(request, align_id):
    """
    Serves query display page
    :param request: HTTP request
    :param align_id: alignment pk
    :return:
    """
    if request.method == 'GET':
        # split sequences in chunks of 80 characters
        length = 80
        align = Alignment.objects.get_alignment(align_id)
        alphabets = {
            "Gapped(ExtendedIUPACProtein(), '-')": 'Protein',
            "Gapped(ExtendedIUPACDNA(), '-')": 'DNA',
        }
        alphabet = alphabets[str(align[0].seq.alphabet)]
        query_seqs = [
            {
                'meta': f.description,
                'seq': [
                    f.seq[i:i + length] for i in range(0, len(f.seq), length)
                    ]
            } for f in align
            ]
        cons_seq = AlignInfo.SummaryInfo(align).gap_consensus()
        consensus = {'meta': 'Consensus', 'seq': [
            cons_seq[i:i + length] for i in range(0, len(cons_seq), length)
            ]}
        return render(
                request,
                'base/query_display.html',
                {'query_seqs': query_seqs, 'seq_type': alphabet, 'consensus': consensus, 'align_id': align_id}
        )

    else:
        return HttpResponseNotAllowed(['GET'])


def align_display(request, align_id):
    """
    serves alignment display page
    :param request: HTTP request
    :param align_id: alignment pk
    :return:
    """
    if request.method == 'GET':
        alignment = Alignment.objects.get_alignment(align_id)
        cons_seq = AlignInfo.SummaryInfo(alignment).gap_consensus()
        cons_seq_annot = {'eq': [0 if c in ['-', 'X'] else 1 for c in cons_seq]}
        cons_seqrec = SeqRecord(
                AlignInfo.SummaryInfo(alignment).gap_consensus(),
                id='consensus',
                name='consensus',
                description='consensus',
                letter_annotations=cons_seq_annot
        )
        alignment.append(cons_seqrec)

        # get longest id for id display width
        id_lengths = [len(a.id) for a in alignment]
        id_lengths.sort(reverse=True)
        id_width = id_lengths[0] + 2

        for a in alignment[:-1]:
            a.letter_annotations = {
                'eq': [1 if (a[i] == cons_seq[i] and cons_seq[i] not in ['-', 'X']) else 0 for i in range(0, len(a))]
            }
        # split sequences in lines of 80 characters
        line_length = 80

        seq_lines = [
            [
                f[i:i + line_length] for f in alignment
                ] for i in range(0, len(alignment[0].seq), line_length)
            ]

        # split lines in blocks of 10 characters
        block_length = 10

        seqs_blocks = [
            [
                [
                    ls.id, [
                        zip(
                            ls.letter_annotations['eq'][j:j + block_length],
                            list(ls[j:j + block_length])
                        ) for j in range(0, len(ls), block_length)
                    ]
                ] for ls in ln] for ln in seq_lines
            ]

        align = {
            'align_seqs': seqs_blocks,
        }
        return render(request, 'base/align_display.html', {'align': align, 'id_width': id_width})
