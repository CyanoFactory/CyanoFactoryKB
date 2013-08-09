'''
CyanoFactory knowledge base tests
'''

from django.test import TestCase, Client

import cyano.models as cmodels

'''
TODO
- import/excel works correctly
'''

class CyanoBaseTest(TestCase):
    def setUp(self):
        self.client = Client()
    
    # via http://djangosnippets.org/snippets/137/
    def GET(self, url, status=200, mimetype="text/html"):
        """Get a URL and require a specific status code before proceeding"""
        response = self.client.get(url)
        self.assertEqual(response.status_code, status, "Status {} (expected: {})".format(response.status_code, status))
        self.assertTrue(response["content-type"].startswith(mimetype), "Mimetype {} (expected: {})".format(response["content-type"], mimetype))
        return response

    def POST(self, url, params, status=200, mimetype="text/html"):
        """Make a POST and require a specific status code before proceeding"""
        response = self.client.post(url, params)
        self.assertEqual(response.status_code, status, "Status {} (expected: {})".format(response.status_code, status))
        self.assertTrue(response["content-type"].startswith(mimetype), "Mimetype {} (expected: {})".format(response["content-type"], mimetype))
        return response

    def assertNotFound(self, url):
        """Fail if found (code != 404)"""
        self.GET(url, status = 404)
    
    def assertForbidden(self, url):
        """Fail if not forbidden (code 403)"""
        self.GET(url, status = 403)

    def assertBadRequest(self, url):
        """Fail if not bad request (code 400)"""
        self.GET(url, status = 400)
    
    def assertRedirect(self, url):
        """Fail if a url does not redirect (code 302)"""
        self.GET(url, status = 302)
    
    def assertRedirectPermanent(self, url):
        """Fail if a url does not permanently redirect (code 301)"""
        self.GET(url, status = 301)
    
    def assertMimeTypeText(self, url, status=200):
        """Fail if mimetype is not text"""
        self.GET(url, status = status, mimetype = "text/plain")
    
    def assertMimeTypeJson(self, url, status=200):
        """Fail if mimetype is not json"""
        self.GET(url, status = status, mimetype = "application/json")
    
    def assertMimeTypeXml(self, url, status=200):
        """Fail if mimetype is not XML"""
        self.GET(url, status = status, mimetype = "application/xml")

class CyanoTest(CyanoBaseTest):
    # Permission list
    fixtures = ['metadata.json']
    
    def setUp(self):
        super(CyanoTest, self).setUp()
        
        # Create users and groups
        from management.commands.autocreateinitial import Command
        Command().handle()

    def test_database_access_works(self):
        cmodels.UserProfile.objects.get(user__username = "admin")
    
    def testTestWorks(self):
        self.assertFalse(True == False, "Not False")
