from django.db import models
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet.IUPAC import ExtendedIUPACProtein, ExtendedIUPACDNA
from Bio.Alphabet import Gapped

# Create your models here.


def save_alignment_to_db(name, data):
    alignment = Alignment.objects.create(name=name)
    for s in data:
        seqrec = Seqrecord.objects.create(
                seq=str(s.seq),
                alphabet=str(s.seq.alphabet),
                seq_id=s.id,
                name=s.name,
                description=s.description
        )
        alignment.seqs.add(seqrec)
    return alignment.pk


def get_multipleseqalignment_object_from_db(pk):
    alignment = Alignment.objects.get(pk=pk).seqs.all()
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
                        # description=a.description
                ) for a in alignment
                ]
    )
    return mul_seq_al


class Seqrecord(models.Model):
    """
    Sequrecords model
    """
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
