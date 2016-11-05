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

from django.test import TestCase
from datetime import datetime, timedelta, timezone
from initialize.init_data import add_old_alignment, set_last_run_at_to_overdue
from base.models import Alignment
from djcelery.models import PeriodicTask
from formalign.settings import CLEAN_OLDER

__author__ = 'Stefan Dieterle'


class InitDataTestCase(TestCase):
    """
    Tests add_to_alignment
    """
    def test_add_to_alignment_adds_the_expected_alignment(self):
        """
        Tests that add_to_alignment adds the expected alignment (slug and
        created)
        """
        add_old_alignment()
        try:
            to_old = Alignment.objects.get(slug='ToOldAlignment01')
            self.assertTrue(to_old.slug == 'ToOldAlignment01')
        except Alignment.DoesNotExist:
            self.fail('Alignment ToOldAlignment01 does not exist')
        self.assertTrue(datetime.now(timezone.utc) - timedelta(days=int(CLEAN_OLDER), hours=1) > to_old.created,
                        'ToOldAlignment created date is not old enough')

    def test_set_last_run_to_overdue_sets_last_run(self):
        """
        Tests that set_last_run_to_overdue sets djcelery_periodictask.last_run
        to a date that is overdue
        """
        PeriodicTask.objects.create(name='remove_old_alignments')
        set_last_run_at_to_overdue()
        clean = PeriodicTask.objects.get(name='remove_old_alignments')
        self.assertTrue(datetime.now(timezone.utc) - timedelta(hours=1, minutes=1) > clean.last_run_at,
                        'remove_old_alignments last_run_at is not old enough')
