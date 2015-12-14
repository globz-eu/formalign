from django.test import TestCase
from django.core.urlresolvers import resolve

from base.views import index, seq_display

__author__ = 'Stefan Dieterle'


class BaseURLsTestCase(TestCase):

    def test_root_url_uses_index_view(self):
        """
        Tests that the root of the site resolves to the correct view function
        :return:
        """
        root = resolve('/')
        self.assertEqual(root.func, index)

    def test_query_seq_url_uses_seq_display_view(self):
        """
        Test that the /query-sequences/ URL resolves to the correct view function
        :return:
        """
        query_seq = resolve('/query-sequences/1')
        self.assertEqual(query_seq.func, seq_display)
