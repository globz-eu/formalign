from django import forms
from helper_funcs.helpers_bio import parse_fasta_alignment
import io

__author__ = 'Stefan Dieterle'

EMPTY_ERROR = 'Please submit an alignment'
FASTA_ERROR = 'Sequence is not FASTA compliant, no ">" as first character'
CHARACTER_ERROR = 'Invalid character in sequence: '
ALIGNMENT_ERROR = 'Alignment invalid, sequences have different lengths'
LESS_THAN_TWO_SEQS_ERROR = 'Submitted data is not a valid alignment, it contains less than 2 sequences'


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

    seq_type_choices = [
        ('Protein', 'Protein'),
        ('DNA', 'DNA')
    ]
    seq_type = forms.ChoiceField(
        widget=forms.RadioSelect(),
        choices=seq_type_choices,
        label='Input sequence type:',
        required=True,
        initial='DNA',
    )

    def clean_align_input(self):
        """
        Returns cleaned and validated alignment sequence data. Validates FASTA for standard FASTA alignment
        (starts with '>', does not contain any invalid characters for protein sequences, all sequences have the same
        length)
        :return: parsed_data = [{'meta': 'sequence meta', 'seq': 'SEQUENCE'} ... ]
        """
        align_input = self.cleaned_data['align_input']
        data = io.StringIO(align_input)

        try:
            align_input = parse_fasta_alignment(data)
        except ValueError:
            raise forms.ValidationError(ALIGNMENT_ERROR)

        if self.cleaned_data['align_input'][0] != '>':
            raise forms.ValidationError(FASTA_ERROR)

        if len(align_input) <= 1:
            raise forms.ValidationError(LESS_THAN_TWO_SEQS_ERROR)

        return align_input

    def clean(self):
        """
        Adds error to align_input field depending on seq_type field
        :return:
        """
        cleaned_data = forms.Form.clean(self)
        seq_type = cleaned_data.get('seq_type')
        align_input = cleaned_data.get('align_input')

        alphabets = {'DNA': set('ACGNT-'), 'Protein': set('ABCDEFGHIKLMNPQRSTUVWXY*-')}
        alphabet = alphabets[seq_type]
        if align_input:
            try:
                for p in align_input:
                    if not alphabet.issuperset(p.seq):
                        self.add_error('align_input', CHARACTER_ERROR + '%s' % p.description)
                        raise ValueError
            except ValueError:
                pass
