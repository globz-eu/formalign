from django.db import models
from randomslugfield import RandomSlugField
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet.IUPAC import ExtendedIUPACProtein, ExtendedIUPACDNA
from Bio.Alphabet import Gapped


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
    created = models.DateTimeField(auto_now_add=True)
    name = models.CharField(max_length=50)
    seqs = models.ManyToManyField(Seqrecord)
    slug = RandomSlugField(length=16)

    objects = AlignmentManager()
