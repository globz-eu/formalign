from django.db import models

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
