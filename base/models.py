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

from django.db import models
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet.IUPAC import ExtendedIUPACProtein, ExtendedIUPACDNA
from Bio.Alphabet import Gapped

__author__ = 'Stefan Dieterle'


class AlignmentManager(models.Manager):
    """
    Manages saving alignments to Alignment model and retrieving data from Alignment model
    """

    def create_alignment(self, name, data):
        """
        Method for saving name to Alignment and saving sequences to Seqrecord model while associating sequences to that
        alignment
        :param name: alignment name
        :param data: list of seqrecords for alignment
        :return: created alignment object
        """
        alignment = self.create(name=name)
        count = 0
        for s in data:
            seqrec = Seqrecord.objects.create(
                display_order=count,
                seq=str(s.seq),
                alphabet=str(s.seq.alphabet),
                seq_id=s.id,
                name=s.name,
                description=s.description
            )
            count += 1
            alignment.seqs.add(seqrec)
        return alignment

    def get_alignment(self, pk):
        """
        Method for managing fetching alignment from db and converting to MultipleSeqAlignment object
        :param pk: alignment pk
        :return: MultipleSeqAlignment object
        """
        alignment = self.get(pk=pk).seqs.all().order_by('display_order')
        alphabets = {
            "Gapped(ExtendedIUPACProtein(), '-')": Gapped(ExtendedIUPACProtein()),
            "Gapped(ExtendedIUPACDNA(), '-')": Gapped(ExtendedIUPACDNA()),
        }
        mul_seq_al = MultipleSeqAlignment(
            [
                SeqRecord(Seq(
                    a.seq, alphabets[a.alphabet]),
                    id=a.seq_id,
                    name=a.name,
                    description=a.description
                ) for a in alignment
                ],
            alphabet=alphabets[alignment[0].alphabet]
        )
        return mul_seq_al


class Seqrecord(models.Model):
    """
    Seqrecords model
    """
    display_order = models.IntegerField(default=0)
    seq = models.TextField()
    alphabet = models.CharField(max_length=50)
    seq_id = models.CharField(max_length=50)
    name = models.CharField(max_length=50)
    description = models.CharField(max_length=200)


class Alignment(models.Model):
    """
    Alignment model
    """
    name = models.CharField(max_length=50)
    seqs = models.ManyToManyField(Seqrecord)

    objects = AlignmentManager()
