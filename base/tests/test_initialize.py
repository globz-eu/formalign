from django.test import TestCase
from datetime import datetime, timedelta, timezone
from initialize.init_data import add_old_alignment, set_last_run_at_to_overdue
from base.models import Alignment
from formalign.settings import CLEAN_OLDER


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
