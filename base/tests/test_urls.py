from django.test import TestCase
from django.core.urlresolvers import resolve

from base.views import index

__author__ = 'Stefan Dieterle'

class BaseURLsTestCase(TestCase):

    def test_root_url_uses_index_view(self):
        """
        Tests that the root of the site resolves to the correct view function
        :return:
        """
        root = resolve('/')
        self.assertEqual(root.func, index)
