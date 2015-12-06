from django import forms
import io
from Bio import SeqIO

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

    def parse_fasta(self, fasta):
        """
        Parses fasta file into list of dicts of metadata and sequences using SeqIO.parse
        :param fasta: fasta string
        :return: fasta_list = [{'meta': 'sequence meta', 'seq': 'SEQUENCE'} ... ]
        """
        fasta_parse = SeqIO.parse(fasta, 'fasta')
        fasta_dict_list = [
            {
               'meta': f.description,
               'seq':  str(f.seq).upper()
            } for f in fasta_parse
            ]
        return fasta_dict_list

    def clean_align_input(self):
        """
        Returns cleaned and validated alignment sequence data. Validates FASTA for standard FASTA alignment
        (starts with '>', does not contain any invalid characters for protein sequences, all sequences have the same
        length)
        :return: parsed_data = [{'meta': 'sequence meta', 'seq': 'SEQUENCE'} ... ]
        """
        if self.cleaned_data['align_input'][0] != '>':
            raise forms.ValidationError(FASTA_ERROR)

        data = io.StringIO(self.cleaned_data['align_input'])
        parsed_data = self.parse_fasta(data)

        if len(parsed_data) <= 1:
            raise forms.ValidationError(LESS_THAN_TWO_SEQS_ERROR)
        alphabet = set('ABCDEFGHIKLMNPQRSTUVWXY*-')
        lengths = []
        for p in parsed_data:
            if not alphabet.issuperset(p['seq']):
                raise forms.ValidationError(CHARACTER_ERROR + '%s' % p['meta'])
            lengths.append(len(p['seq']))

        # check that lengths contains only identical values, count trick is supposed to be faster than using set
        if lengths.count(lengths[0]) != len(lengths):
            raise forms.ValidationError(ALIGNMENT_ERROR)
        return parsed_data
