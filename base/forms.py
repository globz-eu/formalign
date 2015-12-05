from django import forms

__author__ = 'Stefan Dieterle'


EMPTY_ALIGNMENT_SUBMISSION_ERROR = 'Please submit an alignment'


class QueryForm(forms.Form):
    """
    Alignment input form for home page
    """
    align_input = forms.CharField(
        widget=forms.Textarea(
            attrs={
                'placeholder': 'FASTA alignment',
                'class': 'form-control',
            }
        ),
        label='Paste in your alignment in FASTA format:',
        required=True,
        error_messages={'required': 'Please submit an alignment'},
    )
