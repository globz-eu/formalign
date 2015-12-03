from django.test import TestCase, RequestFactory

from base.views import index

__author__ = 'Stefan Dieterle'


class IndexViewTestCase(TestCase):

    def setUp(self):
        self.factory = RequestFactory()

    def test_index_view_basic(self):
        """
        Tests that index view returns a 200 response and uses the correct template
        :return:
        """
        request = self.factory.get('/')
        with self.assertTemplateUsed('base/index.html'):
            response = index(request)
            self.assertEqual(response.status_code, 200)
