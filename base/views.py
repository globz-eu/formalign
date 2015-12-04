from django.shortcuts import render_to_response, redirect, render
from django.http import HttpResponse, HttpResponseRedirect

__author__ = 'Stefan Dieterle'


def index(request):
    return render(request, 'base/index.html')


def seq_display(request):
    if 'align-query' in request.POST and request.POST['align-query']:
        return render(request, 'base/query_display.html', {'query_seqs': request.POST['align-query']})
    # return render(request, 'base/query_display.html', {'query-seqs': request.POST['align-query']})
