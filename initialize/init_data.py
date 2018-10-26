#! python3
import os
import django
from formalign.settings import HEROKU


if HEROKU:
    os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'formalign.settings')
else:
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
    Sets last_run_at in djcelery_periodictask to a date overdue for running
    task
    """
    from datetime import datetime, timedelta, timezone
    from djcelery.models import PeriodicTask

    PeriodicTask.objects.filter(
        name='remove_old_alignments'
    ).update(
        last_run_at=(datetime.now(timezone.utc) - timedelta(hours=1, minutes=1))
    )


if __name__ == '__main__':
    add_old_alignment()
    set_last_run_at_to_overdue()
