from django.shortcuts import render_to_response

__author__ = 'Stefan Dieterle'


def index(request):
    return render_to_response('base/index.html')
