'''
CyanoFactory knowledge base tests
'''

from django.test import TestCase, Client
from django.core import management

import cyano.models as cmodels
from bioparser.genbank import Genbank
from django.core.exceptions import ObjectDoesNotExist

'''
TODO
- import/excel works correctly
'''

SPECIES = "UnitTest-Species"

def setup():
    # Permissions
    management.call_command('loaddata', 'metadata.json')
    # Create users and groups
    management.call_command('autocreateinitial')
    # Create meta data tables
    management.call_command('create_meta')
    # Create simple entry
    g = Genbank(wid = SPECIES,
            user = "management",
            reason = "Unit Testing",
            chromosome = "UnitTest-Chromosome",
            name = "Chromosome for Unit Test")
    with open("../sequences/testcase.gb", "r") as f:
        g.parse(f)
    g.apply()

class CyanoBaseTest(TestCase):
    def setUp(self):
        self.client = Client()
    
    def doLogin(self):
        self.client.login(username = "gabriel", password = "aaa")
        return cmodels.UserProfile.objects.get(user__username = "gabriel")
    
    # via http://djangosnippets.org/snippets/137/
    def GET(self, url, status=200, mimetype="text/html"):
        """Get a URL and require a specific status code before proceeding"""
        if url[0] != "/":
            url = "/" + url
        response = self.client.get(url)
        self.assertEqual(response.status_code, status, "Status {} (expected: {})".format(response.status_code, status))
        self.assertTrue(response["content-type"].startswith(mimetype), "Mimetype {} (expected: {})".format(response["content-type"], mimetype))
        return response

    def POST(self, url, params, status=200, mimetype="text/html"):
        """Make a POST and require a specific status code before proceeding"""
        if url[0] != "/":
            url = "/" + url
        response = self.client.post(url, params)
        self.assertEqual(response.status_code, status, "Status {} (expected: {})".format(response.status_code, status))
        self.assertTrue(response["content-type"].startswith(mimetype), "Mimetype {} (expected: {})".format(response["content-type"], mimetype))
        return response

    def assertOK(self, url):
        """Fail if result code is not 200"""
        return self.GET(url, status = 200)

    def assertNotFound(self, url):
        """Fail if found (code != 404)"""
        return self.GET(url, status = 404)
    
    def assertForbidden(self, url):
        """Fail if not forbidden (code 403)"""
        return self.GET(url, status = 403)

    def assertBadRequest(self, url):
        """Fail if not bad request (code 400)"""
        return self.GET(url, status = 400)
    
    def assertRedirect(self, url, target_url = None, target_status = 200):
        """Fail if a url does not redirect (code 302)"""
        response = self.GET(url, status = 302)
        if target_url:
            if target_url[0] != "/":
                target_url = "/" + target_url
            return self.assertRedirects(response, target_url, 302, target_status)
        return response
    
    def assertRedirectPermanent(self, url, target_url = None, target_status = 200):
        """Fail if a url does not permanently redirect (code 301)"""
        response = self.GET(url, status = 301)
        if target_url:
            if target_url[0] != "/":
                target_url = "/" + target_url
            return self.assertRedirects(response, target_url, 301, target_status)
        return response
    
    def assertMimeTypeText(self, url, status=200):
        """Fail if mimetype is not text"""
        return self.GET(url, status = status, mimetype = "text/plain")
    
    def assertMimeTypeJson(self, url, status=200):
        """Fail if mimetype is not json"""
        return self.GET(url, status = status, mimetype = "application/json")
    
    def assertMimeTypeXml(self, url, status=200):
        """Fail if mimetype is not XML"""
        return self.GET(url, status = status, mimetype = "application/xml")

class CyanoBasicPageTest(CyanoBaseTest):
    def test_append_slash(self):
        """Test if APPEND_SLASH setting is enabled"""
        self.doLogin()
        
        with self.assertTemplateUsed("cyano/login.html"): 
            self.assertRedirectPermanent("login", "login/")
    
    def test_guest(self):
        """Tests if normal page access (as a guest) works"""
        with self.assertTemplateUsed("cyano/index.html"): 
            self.assertOK("/")

        with self.assertTemplateUsed("cyano/logout.html"):
            self.assertOK("logout/")
            self.assertOK(SPECIES + "/logout/")
        
        with self.assertTemplateUsed("cyano/login.html"):
            self.assertOK("login/")
            self.assertOK(SPECIES + "/login/")
        
        #with self.assertTemplateUsed("cyano/about.html"):
        #    self.assertOK("about/")
            
        #with self.assertTemplateUsed("cyano/tutorial.html"):
        #    self.assertOK("tutorial/")

    
    def test_sitemap(self):
        """Tests if sitemaps works"""
        
        with self.assertTemplateUsed("cyano/sitemap.xml"):
            self.assertOK("sitemap.xml")
        
        with self.assertTemplateUsed("cyano/sitemap_toplevel.xml"):
            self.assertOK("sitemap_toplevel.xml")
        
        with self.assertTemplateUsed("cyano/sitemap_species.xml"):
            self.assertOK("sitemap/%s.xml" % (SPECIES))
        
        with self.assertTemplateUsed("cyano/error.html"):
            self.assertNotFound("sitemap/Wrong-Species.xml")

class CyanoTest(CyanoBaseTest):
    def test_logic_works(self):
        self.assertFalse(True == False, "Hardware error")
    
    def test_permissions_species(self):
        self.assertForbidden(SPECIES + "/")
        
        user = self.doLogin()

        self.assertForbidden(SPECIES + "/")

        entry = cmodels.Species.objects.for_wid(SPECIES)
        user.add_allow_permission(entry, "READ_NORMAL")

        with self.assertTemplateUsed("cyano/species.html"): 
            self.assertOK(SPECIES + "/")
        
        user.add_deny_permission(entry, "READ_NORMAL")
        
        self.assertForbidden(SPECIES + "/")
    
    def test_empty_permission_deleted(self):
        user = self.doLogin()
        entry = cmodels.Species.objects.for_wid(SPECIES)
        
        with self.assertRaises(ObjectDoesNotExist):
            cmodels.UserPermission.objects.get(user = user, entry = entry)
        
        user.add_allow_permission(entry, "READ_NORMAL")
        cmodels.UserPermission.objects.get(user = user, entry = entry)
        user.delete_allow_permission(entry, "READ_NORMAL")
        
        with self.assertRaises(ObjectDoesNotExist):
            cmodels.UserPermission.objects.get(user = user, entry = entry)
    
    def test_species_not_found(self):
        with self.assertTemplateUsed("cyano/error.html"):
            self.assertNotFound("/Wrong-Species/")

