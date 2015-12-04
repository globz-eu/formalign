from django.test import TestCase

__author__ = 'Stefan Dieterle'


class IndexViewTestCase(TestCase):

    def test_index_view_basic(self):
        """
        Tests that index view returns a 200 response and uses the correct template
        :return:
        """
        response = self.client.get('/')
        self.assertTemplateUsed('base/index.html')
        self.assertEqual(response.status_code, 200)


class SeqDisplayTestCase(TestCase):

    def test_display_page_uses_display_seq_template(self):
        """
        Tests that seq_display view returns a 200 response on a POST request and uses the correct template
        :return:
        """
        response = self.client.post('/query-sequences/', {'align-query': '>New_seq\nGLOK'})
        self.assertTemplateUsed('base/query_display.html')
        self.assertEqual(response.status_code, 200)

    def test_display_page_displays_query_seq(self):
        """
        Tests that seq_display displays the query on a valid POST request
        :return:
        """
        response = self.client.post('/query-sequences/', {'align-query': '>New_seq\nGLOK'})
        self.assertIn('&gt;New_seq\nGLOK', response.content.decode('utf-8'), 1)
