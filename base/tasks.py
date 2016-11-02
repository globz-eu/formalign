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

from __future__ import absolute_import
from celery import shared_task
from datetime import datetime, timedelta, timezone
from formalign.settings import CLEAN_OLDER

from base.models import Alignment, Seqrecord

__author__ = 'Stefan Dieterle'


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
