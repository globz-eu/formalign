"""
base app views
"""

import io

from Bio.Alphabet import Gapped
from Bio.Alphabet.IUPAC import ExtendedIUPACProtein
from django.http import HttpResponseNotAllowed, Http404
from django.shortcuts import render, redirect

from base.forms import QueryForm
from base.models import Alignment
from helper_funcs.bio.helpers import consensus_add, parse_fasta_alignment
from helper_funcs.helpers_format import split_lines, annotate
from helper_funcs.helpers_test import file_to_string


def index(request):
    """
    Serves home page
    :param request: HTTP request
    :return: HttpResponse object
    """
    if request.method == 'GET':
        form = QueryForm()
        return render(request, 'base/index.html', {'form': form})
    if request.method == 'POST':
        if request.POST['custom_data'] == 'custom':
            form = QueryForm(request.POST)
            if form.is_valid():
                align = form.cleaned_data['align_input']
                save_align = Alignment.objects.create_alignment('name', align)
                slug = save_align.slug
                return redirect('/query-sequences/' + str(slug) + '/')
            return render(request, 'base/index.html', {'form': form})
        if request.POST['custom_data'] == 'demo':
            align_input = io.StringIO(file_to_string('ser_thr_kinase_family.fasta'))
            align = parse_fasta_alignment(align_input)
            for sequence in align:
                sequence.seq.alphabet = Gapped(ExtendedIUPACProtein())
            save_align = Alignment.objects.create_alignment('ser_thr_kinase_family', align)
            slug = save_align.slug
            return redirect('/query-sequences/' + str(slug) + '/')
        form = QueryForm()
        return render(request, 'base/index.html', {'form': form})
    return HttpResponseNotAllowed(['POST', 'GET'])


def seq_display(request, align_slug):
    """
    Serves query display page
    :param request: HTTP request
    :param align_id: alignment pk
    :return:HttpResponse or HttpResponseNotAllowed object
    """
    if request.method == 'GET':
        # split sequences in chunks of 80 characters
        try:
            align_id = Alignment.objects.get(slug=align_slug).pk
        except Alignment.DoesNotExist:  # pylint: disable=E1101
            raise Http404('Alignment does not exist')
        alignment = consensus_add(Alignment.objects.get_alignment(align_id))
        alphabets = {
            "Gapped(ExtendedIUPACProtein(), '-')": 'Protein',
            "Gapped(ExtendedIUPACDNA(), '-')": 'DNA',
        }
        alphabet = alphabets[str(alignment[0].seq.alphabet)]

        query_seqs = split_lines(alignment, line_length=80, split_type='sequence')

        return render(
            request,
            'base/query_display.html',
            {'query_seqs': query_seqs, 'seq_type': alphabet, 'align_id': align_slug}
        )

    return HttpResponseNotAllowed(['GET'])


def align_display(request, align_slug):
    """
    serves alignment display page
    :param request: HTTP request
    :param align_slug: alignment slug value
    :return: HttpResponse object
    """
    if request.method == 'GET':
        try:
            align_id = Alignment.objects.get(slug=align_slug).pk
        except Alignment.DoesNotExist:  # pylint: disable=E1101
            raise Http404('Alignment does not exist')

        alignment_fetch = Alignment.objects.get_alignment(align_id)
        alignment = consensus_add(alignment_fetch)

        alignment = annotate(alignment)

        # get longest id for id display width
        id_lengths = [len(a.id) for a in alignment]
        id_lengths.sort(reverse=True)
        id_width = id_lengths[0] + 2

        # split sequences in lines of 80 characters
        seq_lines = split_lines(alignment, line_length=80, split_type='alignment')

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

        return render(
            request, 'base/align_display.html',
            {'align': align, 'id_width': id_width * .6, 'total_width': (id_width + 80) * .5}
        )

    return HttpResponseNotAllowed(['GET'])
