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

from behave import given, when
from formalign.settings import SERVER_URL
from initialize.init_data import add_old_alignment
from base.tasks import clean_alignments
from django.test.utils import override_settings

__author__ = 'Stefan Dieterle'


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


@when(r'a user visits the URL "(?P<url>[^"]*)"')
def visit_url(context, url):
    """
    visits the given URL from the base server URL
    :param context: behave context
    :param url: URL to visit
    """
    if SERVER_URL == 'liveserver':
        context.home_url = context.base_url
    else:
        context.home_url = SERVER_URL
    context.r = context.client.get(context.home_url + url)
