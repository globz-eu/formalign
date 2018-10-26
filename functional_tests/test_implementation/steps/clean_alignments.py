from behave import given, use_step_matcher
from formalign.settings import SERVER_URL
from initialize.init_data import add_old_alignment
from base.tasks import clean_alignments
from django.test.utils import override_settings


use_step_matcher('re')


@given(r'an old alignment with slug "ToOldAlignment01" was present')
def initialize_with_old_alignment(context):
    """
    creates an old alignment
    :param context: behave context
    """
    if SERVER_URL == 'liveserver':
        add_old_alignment()


@given(r'clean_alignments task has been run')
def run_clean_alignments(context):
    """
    runs clean_alignments task
    :param context: behave context
    """
    if SERVER_URL == 'liveserver':
        @override_settings(CELERY_EAGER_PROPAGATES_EXCEPTIONS=True, CELERY_ALWAYS_EAGER=True)
        def clean():
            clean_alignments()
    else:
        def clean():
            pass
    clean()
