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
                {'query_seqs': form.cleaned_data['align_input']}
            )
        else:
            return render(request, 'base/index.html', {'form': form})
