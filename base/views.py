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

import io

from django.shortcuts import render, redirect
from django.http import HttpResponseNotAllowed

from Bio.Alphabet.IUPAC import ExtendedIUPACProtein
from Bio.Alphabet import Gapped

from base.models import Alignment
from base.forms import QueryForm

from helper_funcs.helpers_bio import consensus_add, parse_fasta_alignment
from helper_funcs.helpers_format import split_lines, annotate
from helper_funcs.helpers_test import file_to_string

__author__ = 'Stefan Dieterle'


def index(request):
    """
    Serves home page
    :param request: HTTP request
    :return: HttpResponse object
    """
    if request.method == 'GET':
        form = QueryForm()
        return render(request, 'base/index.html', {'form': form})
    elif request.method == 'POST':
        if request.POST['custom_data'] == 'custom':
            form = QueryForm(request.POST)
            if form.is_valid():
                align = form.cleaned_data['align_input']
                save_align = Alignment.objects.create_alignment('name', align)
                slug = save_align.slug
                return redirect('/query-sequences/' + str(slug) + '/')
            else:
                return render(request, 'base/index.html', {'form': form})
        elif request.POST['custom_data'] == 'demo':
            align_input = io.StringIO(file_to_string('ser_thr_kinase_family.fasta'))
            align = parse_fasta_alignment(align_input)
            for d in align:
                d.seq.alphabet = Gapped(ExtendedIUPACProtein())
            save_align = Alignment.objects.create_alignment('ser_thr_kinase_family', align)
            slug = save_align.slug
            return redirect('/query-sequences/' + str(slug) + '/')
        else:
            form = QueryForm()
            return render(request, 'base/index.html', {'form': form})
    else:
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
        align_id = Alignment.objects.get(slug=align_slug).pk
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

    else:
        return HttpResponseNotAllowed(['GET'])


def align_display(request, align_slug):
    """
    serves alignment display page
    :param request: HTTP request
    :param align_id: alignment pk
    :return: HttpResponse object
    """
    if request.method == 'GET':
        align_id = Alignment.objects.get(slug=align_slug).pk
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

        r = render(
            request, 'base/align_display.html',
            {'align': align, 'id_width': id_width * .6, 'total_width': (id_width + 80) * .5}
        )
        return r
