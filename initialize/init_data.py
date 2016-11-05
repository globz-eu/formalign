#! python3
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

import os
import django

__author__ = 'Stefan Dieterle'

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'settings_admin')
django.setup()


def add_old_alignment():
    """
    Adds an alignment with a creation date older than CLEAN_OLDER days
    """
    import io
    import sys
    from datetime import datetime, timedelta, timezone
    from helper_funcs.bio.helpers import parse_fasta_alignment
    from helper_funcs.helpers_test import file_to_string
    from base.models import Alignment
    from Bio.Alphabet.IUPAC import ExtendedIUPACProtein
    from Bio.Alphabet import Gapped
    from formalign.settings import CLEAN_OLDER

    name = 'A. tha. SPA family alignment'
    align_input = io.StringIO(file_to_string('spa_protein_alignment.fasta'))
    data = parse_fasta_alignment(align_input)
    alphabet = Gapped(ExtendedIUPACProtein())
    for a in data:
        a.seq.alphabet = alphabet
    data._alphabet = alphabet
    try:
        Alignment.objects.get(slug='ToOldAlignment01')
        Alignment.objects.get(slug='ToOldAlignment01').delete()
        print('removed alignment', file=sys.stderr)
    except Alignment.DoesNotExist:
        pass

    alignment = Alignment.objects.create_alignment(name, data)
    Alignment.objects.filter(
        slug=alignment.slug
    ).update(
        created=(datetime.now(timezone.utc) - timedelta(days=int(CLEAN_OLDER), hours=1)),
        slug='ToOldAlignment01'
    )


def set_last_run_at_to_overdue():
    """
    Sets last_run_at in djcelery_periodictask to a date overdue for running task
    """
    from datetime import datetime, timedelta, timezone
    from djcelery.models import PeriodicTask

    PeriodicTask.objects.filter(
        name='remove_more_than_week_old_alignments'
    ).update(
        last_run_at=(datetime.now(timezone.utc) - timedelta(hours=1, minutes=1))
    )


if __name__ == '__main__':
    add_old_alignment()
    set_last_run_at_to_overdue()
