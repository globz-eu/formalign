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

from Bio.Alphabet import Gapped
from Bio.Alphabet.IUPAC import ExtendedIUPACProtein, ExtendedIUPACDNA
from django import forms
from django.utils.safestring import mark_safe
from helper_funcs.bio.helpers import parse_fasta_alignment

__author__ = 'Stefan Dieterle'

EMPTY_ERROR = 'Please submit an alignment'
FORMAT_ERROR = """The server could not figure out what format this is,
please double check your input or try a different format"""
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
                'placeholder': 'Alignment (FASTA, clustalw, stockholm or phylip)',
                'class': 'form-control',
            }
        ),
        label=mark_safe('Paste in your alignment:<br>(FASTA, clustalw, stockholm or phylip)'),
        label_suffix='',
        required=True,
        error_messages={'required': 'Please submit an alignment'},
    )

    seq_type_choices = [
        ('protein', 'Protein'),
        ('DNA', 'DNA')
    ]
    seq_type = forms.ChoiceField(
        widget=forms.RadioSelect(),
        choices=seq_type_choices,
        label='Input sequence type:',
        required=True,
        initial='DNA',
    )

    cons_type_choices = [
        ('identity', 'Identity'),
        ('subs_matrix', 'Substitution Matrix')
    ]
    cons_type = forms.ChoiceField(
        widget=forms.RadioSelect(),
        choices=cons_type_choices,
        label='Consensus type:',
        required=True,
        initial='identity',
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

        if self.cleaned_data['align_input'][0] != '>':
            raise forms.ValidationError(FORMAT_ERROR)

        try:
            align_input = parse_fasta_alignment(data)
        except ValueError:
            raise forms.ValidationError(ALIGNMENT_ERROR)

        if len(align_input) <= 1:
            raise forms.ValidationError(LESS_THAN_TWO_SEQS_ERROR)

        return align_input

    def clean(self):
        """
        Adds error to align_input field depending on seq_type field, adds alphabets to sequences and alignment
        :return:
        """
        cleaned_data = forms.Form.clean(self)
        seq_type = cleaned_data.get('seq_type')
        alphabets = {'DNA': Gapped(ExtendedIUPACDNA()), 'protein': Gapped(ExtendedIUPACProtein())}

        align_input = cleaned_data.get('align_input')
        if align_input:
            for a in cleaned_data['align_input']:
                a.seq.alphabet = alphabets[seq_type]
            cleaned_data['align_input']._alphabet = alphabets[seq_type]
        alphabet = set(alphabets[seq_type].letters)
        if align_input:
            try:
                for p in align_input:
                    if not alphabet.issuperset(p.seq):

                        raise CharacterError(p.description)
            except CharacterError as v_err:
                self.add_error('align_input', v_err.message + '%s' % v_err.seq)


class CharacterError(Exception):
    """
    Custom error class for invalid character in sequence
    """
    def __init__(self, seq):
        self.seq = seq
        self.message = CHARACTER_ERROR
