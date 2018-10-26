from __future__ import absolute_import
from celery import shared_task
from datetime import datetime, timedelta, timezone
from formalign.settings import CLEAN_OLDER

from base.models import Alignment, Seqrecord


@shared_task
def clean_alignments():
    """
    Removes alignments older than CLEAN_OLDER (in days)
    :return: list of names of the removed alignments
    """
    days = timedelta(days=int(CLEAN_OLDER))
    old_alignments = Alignment.objects.filter(created__lte=(datetime.now(timezone.utc) - days))
    old_alignments_names = []
    old_seq_ids = []
    for old in old_alignments:
        old_seq_ids.extend([o.pk for o in old.seqs.all()])
        old_alignments_names.append(old.name)
    Seqrecord.objects.filter(pk__in=old_seq_ids).delete()
    old_alignments.delete()
    return old_alignments_names
