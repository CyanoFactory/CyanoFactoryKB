'''
CyanoFactory knowledge base tests
'''

from django.test import TestCase, Client
from django.core import management

import cyano.models as cmodels
from bioparser.genbank import Genbank
from django.core.exceptions import ObjectDoesNotExist

import settings

'''
TODO
- import/excel works correctly
'''

SPECIES = "UnitTest-Species"
SPECIES2 = "Other-Species"
MODEL = "Gene"
ITEM = "sll7106"

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
    
    species = cmodels.Species(wid = "Other-Species", name = "Unit Test species 2")
    rev = cmodels.RevisionDetail(
        user = cmodels.UserProfile.objects.get(user__username = "management"),
        reason = "Unit Testing")
    rev.save()
    species.save(rev)

class CyanoBaseTest(TestCase):
    def setUp(self):
        self.client = Client()
    
    def doLogin(self):
        self.client.login(username = "gabriel", password = "aaa")
        return cmodels.UserProfile.objects.get(user__username = "gabriel")
    
    def doAdminLogin(self):
        self.client.login(username = "admin", password = "aaa")
        return cmodels.UserProfile.objects.get(user__username = "admin")

    def doGuestLogin(self):
        """Fake login returning the guest handle"""
        return cmodels.UserProfile.objects.get(user__username = "guest")
    
    def getSpecies(self, second = False):
        return cmodels.Species.objects.for_wid(SPECIES2 if second else SPECIES)
    
    def getItem(self):
        return getattr(cmodels, MODEL).objects.for_species(SPECIES).for_wid(ITEM)

    def checkMany(self, li, assert_func):
        for l in li:
            assert_func(l)

    # via http://djangosnippets.org/snippets/137/
    def GET(self, url, params = {}, status=200, mimetype="text/html", follow = False):
        """Get a URL and require a specific status code before proceeding"""
        if url[0] != "/":
            url = "/" + url
        response = self.client.get(url, params, follow = follow)

        if response.status_code != status:
            # print note that login was required when a test fails
            if "location" in response._headers:
                if "login" in response._headers["location"][1]:
                    print "Login required."

        self.assertEqual(response.status_code, status, "\nURL: {}\nStatus {} (expected: {})".format(url, response.status_code, status))
            
        self.assertTrue(response["content-type"].startswith(mimetype), "Mimetype {} (expected: {})".format(response["content-type"], mimetype))     
        return response

    def POST(self, url, params = {}, status=200, mimetype="text/html"):
        """Make a POST and require a specific status code before proceeding"""
        if url[0] != "/":
            url = "/" + url
        response = self.client.post(url, params)
        self.assertEqual(response.status_code, status, "\nURL: {}\nStatus {} (expected: {})".format(url, response.status_code, status))
        self.assertTrue(response["content-type"].startswith(mimetype), "Mimetype {} (expected: {})".format(response["content-type"], mimetype))
        return response

    def assertOK(self, url, data = {}, follow = False):
        """Fail if result code is not 200"""
        return self.GET(url, params = data, status = 200, follow = follow)
    
    def assertPOSTOK(self, url, data = {}, follow = False):
        """Fail if result code is not 200"""
        return self.POST(url, params = data, status = 200)

    def assertNotFound(self, url, follow = False):
        """Fail if found (code != 404)"""
        return self.GET(url, status = 404, follow = follow)
    
    def assertForbidden(self, url, follow = False):
        """Fail if not forbidden (code 403)"""
        return self.GET(url, status = 403, follow = follow)

    def assertBadRequest(self, url, follow = False):
        """Fail if not bad request (code 400)"""
        return self.GET(url, status = 400, follow = follow)
    
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
    
    def assertMimeTypeText(self, url, status=200, follow = False):
        """Fail if mimetype is not text"""
        return self.GET(url, status = status, mimetype = "text/plain", follow = follow)
    
    def assertMimeTypeJson(self, url, status=200, follow = False):
        """Fail if mimetype is not json"""
        return self.GET(url, status = status, mimetype = "application/json", follow = follow)
    
    def assertMimeTypeXml(self, url, status=200, follow = False):
        """Fail if mimetype is not XML"""
        return self.GET(url, status = status, mimetype = "application/xml", follow = follow)
    
    def assertLoginRequired(self, url):
        """Tests if a view redirects to login page"""
        # Not for testing the decorator, only if the decorator is set
        
        with self.assertTemplateUsed("cyano/login.html"):
            self.assertOK(url, follow = True)

class CyanoBasicTest(CyanoBaseTest):
    def test_logic_works(self):
        self.assertFalse(True == False, "Hardware error")

    def test_append_slash(self):
        """Test if APPEND_SLASH setting is enabled"""
        self.doLogin()
        
        with self.assertTemplateUsed("cyano/login.html"): 
            self.assertRedirectPermanent("login", "login/")

class CyanoGuestUserTestBase(CyanoBaseTest):
    def _test_index(self):
        with self.assertTemplateUsed("cyano/index.html"): 
            self.assertOK("/")

    def _test_loginout(self, s):        
        with self.assertTemplateUsed("cyano/%s.html" % (s)):
            self.assertOK("%s/" % (s))
            self.assertOK(SPECIES + "/%s/" % (s))

        self.assertNotFound("%s/invalid/" % (s))
        self.assertNotFound(SPECIES + "/%s/invalid/" % (s))
        self.assertNotFound(SPECIES + "/model/%s/" % (s))

    def _test_search(self, google):
        test_urls = ["%ssearch/",
                     "%ssearch/?engine=haystack",
                     "%ssearch/?engine=",
                     "%ssearch/?engine=abcdef",
                     "%ssearch/?q=query&engine=haystack",
                     "%ssearch/?q=query"]
    
        test_google_urls = ["%ssearch/?engine=google",
                            "%ssearch/?engine=google&query=test"]
        
        with self.assertTemplateUsed("cyano/search.html"):
            for url in test_urls:
                self.assertOK(url % (""))
                self.assertOK(url % (SPECIES + "/"))
        
        with self.assertTemplateUsed("cyano/googleSearch.html" if google else "cyano/search.html"):
            for url in test_google_urls:
                self.assertOK(url % (""))
                self.assertOK(url % (SPECIES + "/"))
                
    def _test_sitemap(self):
        self.assertNotFound("sitemap.xml/")
        self.assertNotFound("sitemap_toplevel.xml/")
        self.assertNotFound("sitemap/%s.xml/" % (SPECIES))
        
        with self.assertTemplateUsed("cyano/sitemap.xml"):
            self.assertOK("sitemap.xml")
        
        with self.assertTemplateUsed("cyano/sitemap_toplevel.xml"):
            self.assertOK("sitemap_toplevel.xml")
        
        with self.assertTemplateUsed("cyano/sitemap_species.xml"):
            self.assertOK("sitemap/%s.xml" % (SPECIES))
        
        with self.assertTemplateUsed("cyano/error.html"):
            self.assertNotFound("sitemap/Wrong-Species.xml")

class CyanoGuestTest(CyanoGuestUserTestBase):
    def test_guest_index(self):
        """Visits the main page as guest"""
        self._test_index()

    def test_guest_login(self):
        """Visit login page as guest"""
        self._test_loginout("login")
    
    def test_guest_logout(self):
        """Visit logout page as guest"""
        self._test_loginout("logout")
    
    def _test_basic_uri_stuff(self, s):
        self.assertLoginRequired("%s/" % (s))
        self.assertLoginRequired(SPECIES + "/%s/" % (s))
        
        self.assertNotFound("%s/invalid/" % (s))
        self.assertNotFound(SPECIES + "/%s/invalid/" % (s))
        self.assertNotFound(SPECIES + "/%s/users/" % (s))

    def test_guest_users(self):
        """Visit users page as guest"""
        self._test_basic_uri_stuff("users")
    
    def test_guest_user(self):
        """Visit user profile pages as guest"""
        self._test_basic_uri_stuff("user/admin")
        self._test_basic_uri_stuff("user/invalid")
    
    def test_guest_job(self):
        """Visit jobs page as guest"""
        self._test_basic_uri_stuff("jobs")
    
    def test_guest_search_google_on(self):
        """Tests search functions (google search enabled) as guest"""
        settings.GOOGLE_SEARCH_ENABLED = True

        self._test_search(google = True)

    def test_guest_search_google_off(self):
        """Tests search functions (google search disabled) as guest"""
        settings.GOOGLE_SEARCH_ENABLED = False

        self._test_search(google = False)

    def test_guest_species(self):
        """Visit species page as guest"""
        self.assertForbidden(SPECIES + "/")
        self.doGuestLogin().add_allow_permission(self.getSpecies(), "READ_NORMAL")
        
        with self.assertTemplateUsed("cyano/species.html"):
            self.assertOK(SPECIES + "/")
    
    def test_guest_listing(self):
        """Visit listing page as guest"""
        url = SPECIES + "/" + MODEL + "/"
        
        self.assertForbidden(url)
        self.assertNotFound(SPECIES + "/" + "Invalid" + "/")
        
        self.doGuestLogin().add_allow_permission(self.getSpecies(), "READ_NORMAL")
        
        self.assertNotFound(SPECIES + "/" + "Invalid" + "/")
        
        with self.assertTemplateUsed("cyano/list.html"):
            self.assertOK(url)
        
        self.assertNotFound(SPECIES + "/Species/")
        self.assertNotFound(SPECIES + "/Entry/")

    def test_guest_detail(self):
        """Visit detail page as guest"""
        url = SPECIES + "/" + MODEL + "/" + ITEM + "/"
        url2 = SPECIES2 + "/" + MODEL + "/" + ITEM + "/"
        self.assertForbidden(url)
        self.assertNotFound(url2)

        self.assertNotFound(SPECIES + "/" + "Invalid" + "/" + ITEM + "/")
        self.assertNotFound(SPECIES + "/" + MODEL + "/" + "Invalid" + "/")

        self.doGuestLogin().add_allow_permission(self.getSpecies(), "READ_NORMAL")
        self.doGuestLogin().add_allow_permission(self.getSpecies(second = True), "READ_NORMAL")
        
        self.assertNotFound(url2)
        self.assertNotFound(SPECIES + "/" + "Invalid" + "/" + ITEM + "/")
        self.assertNotFound(SPECIES + "/" + MODEL + "/" + "Invalid" + "/")
        
        with self.assertTemplateUsed("cyano/detail.html"):
            self.assertOK(url)

    def test_guest_add(self):
        """Visit add page as guest"""
        self.assertNotFound("Species/")
        
        self.assertLoginRequired("Species/add/")
        
        self.assertLoginRequired(SPECIES + "/" + MODEL + "/add/")
        
        self.assertNotFound(SPECIES + "/add/")
        self.assertNotFound(SPECIES + "/" + MODEL + "/" + ITEM + "/add/")
        self.assertNotFound(SPECIES + "/" + MODEL + "/" + ITEM + "/add/invalid/")
        self.assertNotFound(SPECIES + "/add/invalid/")
        self.assertNotFound(SPECIES + "/" + MODEL + "/add/invalid/")
    
    def test_guest_add_permission(self):
        """Visit add page as guest with permission"""
        self.doGuestLogin().add_allow_permission(self.getSpecies(), "WRITE_NORMAL")

        self.assertLoginRequired(SPECIES + "/" + MODEL + "/add/")
    
    def _test_editdelete_guest(self, s):
        self.assertNotFound("%s/" % (s))
        self.assertLoginRequired("Species/%s/" % (s))

        self.assertNotFound(SPECIES + "/" + MODEL + "/%s/" % (s))
        self.assertNotFound(SPECIES + "/" + MODEL + "/" + ITEM + "/%s/invalid/" % (s))
        self.assertNotFound(SPECIES + "/%s/invalid/" % (s))
        self.assertNotFound(SPECIES + "/" + MODEL + "/%s/invalid/" % (s))
        
        self.assertLoginRequired(SPECIES + "/%s/" % (s))
        self.assertLoginRequired(SPECIES + "/" + MODEL + "/" + ITEM + "/%s/" % (s))
    
    def test_guest_edit(self):
        """Visit edit page as guest"""
        self._test_editdelete_guest("edit")
    
    def test_guest_edit_permission(self):
        """Visit edit page as guest with permission"""
        self.doGuestLogin().add_allow_permission(self.getSpecies(), "WRITE_NORMAL")

        self.assertLoginRequired(SPECIES + "/edit/")
    
    def test_guest_delete(self):
        """Visit delete page as guest"""
        self._test_editdelete_guest("delete")
    
    def test_guest_delete_permission(self):
        """Visit delete page as guest with permission"""
        self.doGuestLogin().add_allow_permission(self.getSpecies(), "WRITE_DELETE")

        self.assertLoginRequired(SPECIES + "/delete/")

    def test_guest_exportData(self):
        """Visit export page as guest"""
        self.assertForbidden(SPECIES + "/export/")
        
        self.doGuestLogin().add_allow_permission(self.getSpecies(), "READ_NORMAL")
        
        with self.assertTemplateUsed("cyano/exportDataForm.html"):
            self.assertOK(SPECIES + "/export/")
        
        self.assertNotFound("export/")
        self.assertNotFound("invalid/export/")
        self.assertNotFound(SPECIES + "/model/export/")

    def test_guest_import(self):
        """Visit import page as guest"""
        self.assertNotFound("import/")
        
        self.assertLoginRequired("/import/data/")
        self.assertLoginRequired("/import/species/")
        
        self.assertLoginRequired(SPECIES + "/import/data/")
        self.assertLoginRequired(SPECIES + "/import/species/")
        
        self.assertNotFound("/import/data/invalid/")
        self.assertNotFound("/import/species/invalid/")
        self.assertNotFound(SPECIES + "/import/data/invalid/")
        self.assertNotFound(SPECIES + "/import/species/invalid/")

    def test_guest_history(self):
        """Visit history pages as guest"""
        self.assertNotFound("history/")
        
        self.assertForbidden(SPECIES + "/history/")
        self.assertForbidden(SPECIES + "/" + MODEL + "/history/")
        self.assertForbidden(SPECIES + "/" + MODEL + "/" + ITEM + "/history/")
        
        self.doGuestLogin().add_allow_permission(self.getSpecies(), "READ_HISTORY")
        
        with self.assertTemplateUsed("cyano/history.html"):
            self.assertOK(SPECIES + "/history/")
            self.assertOK(SPECIES + "/" + MODEL + "/history/")
            self.assertOK(SPECIES + "/" + MODEL + "/" + ITEM + "/history/")
    
    def test_guest_history_detail(self):
        """Visit history detail page as guest"""
        self.assertNotFound(SPECIES + "/" + MODEL + "/" + ITEM + "/history/a/")
        self.assertNotFound(SPECIES + "/" + MODEL + "/" + ITEM + "/history/1/a/")
        
        self.assertForbidden(SPECIES + "/" + MODEL + "/" + ITEM + "/history/1/")
        
        self.doGuestLogin().add_allow_permission(self.getSpecies(), "READ_HISTORY")
        
        with self.assertTemplateUsed("cyano/history_detail.html"):
            self.assertOK(SPECIES + "/" + MODEL + "/" + ITEM + "/history/1/")

    def test_guest_permission(self):
        """Visit permission view page as guest"""
        self.assertNotFound("permission/")
        
        self.assertForbidden(SPECIES + "/permission/")
        self.assertNotFound(SPECIES + "/" + MODEL + "/permission/")
        self.assertForbidden(SPECIES + "/" + MODEL + "/" + ITEM + "/permission/")
        
        self.doGuestLogin().add_allow_permission(self.getSpecies(), "READ_PERMISSION")
        
        with self.assertTemplateUsed("cyano/permission.html"):
            self.assertOK(SPECIES + "/permission/")
            self.assertOK(SPECIES + "/" + MODEL + "/" + ITEM + "/permission/")

    def test_guest_permission_edit(self):
        """Visit permission edit page as guest"""
        # Special case: Matches species_wid/edit
        self.assertLoginRequired("permission/edit/")
        # Special case: Matches species_wid/model/edit
        self.assertLoginRequired(SPECIES + "/" + MODEL + "/permission/edit/")
        
        self.assertForbidden(SPECIES + "/permission/edit/")
        self.assertForbidden(SPECIES + "/" + MODEL + "/" + ITEM + "/permission/edit/")
        
        self.doGuestLogin().add_allow_permission(self.getSpecies(), "WRITE_PERMISSION")
        
        self.assertForbidden(SPECIES + "/permission/edit/")
        self.assertForbidden(SPECIES + "/" + MODEL + "/" + ITEM + "/permission/edit/")
        
        self.doGuestLogin().add_allow_permission(self.getSpecies(), "READ_PERMISSION")
        
        with self.assertTemplateUsed("cyano/permission_edit.html"):
            self.assertOK(SPECIES + "/permission/edit/")
            self.assertOK(SPECIES + "/" + MODEL + "/" + ITEM + "/permission/edit/")

    def test_guest_sitemap(self):
        """Visit sitemaps as guest"""
        self._test_sitemap()

class CyanoUserTest(CyanoGuestUserTestBase):
    def setUp(self):
        super(CyanoUserTest, self).setUp()
        self.doLogin()
    
    def test_user_index(self):
        """Visits the main page as logged in user"""
        self._test_index()

    def test_user_login(self):
        """Visit login page as logged in user"""
        self._test_loginout("login")
    
    def test_user_logout(self):
        """Visit logout page as logged in user"""
        self._test_loginout("logout")

    def test_user_users(self):
        """Visit users page as logged in user"""
        with self.assertTemplateUsed("cyano/users.html"):
            self.assertOK("users/")
            self.assertOK(SPECIES + "/users/")
    
    def test_user_user(self):
        """Visit user profile pages as logged in user"""
        with self.assertTemplateUsed("cyano/user.html"):
            self.assertOK("user/admin/")
            self.assertOK(SPECIES + "/user/admin/")
            
        self.assertNotFound("user/invalid/")
        self.assertNotFound(SPECIES + "/user/invalid/")
    
    def test_user_job(self):
        """Visit jobs page as logged in user"""
        with self.assertTemplateUsed("cyano/jobs.html"):
            self.assertOK("jobs/")

    def test_user_search_google_on(self):
        """Tests search functions (google search enabled)"""
        settings.GOOGLE_SEARCH_ENABLED = True
        
        self._test_search(google = True)

    def test_user_search_google_off(self):
        """Tests search functions (google search disabled)"""
        settings.GOOGLE_SEARCH_ENABLED = False
        
        self._test_search(google = False)

    def test_user_species(self):
        """Visit species page as logged in user"""
        self.assertForbidden(SPECIES + "/")
        self.doLogin().add_allow_permission(self.getSpecies(), "READ_NORMAL")
        
        with self.assertTemplateUsed("cyano/species.html"):
            self.assertOK(SPECIES + "/")
    
    def test_user_listing(self):
        """Visit listing page as logged in user"""
        url = SPECIES + "/" + MODEL + "/"
        
        self.assertForbidden(url)

        self.doLogin().add_allow_permission(self.getSpecies(), "READ_NORMAL")
        
        with self.assertTemplateUsed("cyano/list.html"):
            self.assertOK(url)

    def test_user_detail(self):
        """Visit detail page as logged in user"""
        url = SPECIES + "/" + MODEL + "/" + ITEM + "/"
        self.assertForbidden(url)

        self.doLogin().add_allow_permission(self.getSpecies(), "READ_NORMAL")
        self.doLogin().add_allow_permission(self.getSpecies(second = True), "READ_NORMAL")

        with self.assertTemplateUsed("cyano/detail.html"):
            self.assertOK(url)
            self.assertOK(SPECIES + "/Chromosome/UnitTest-Chromosome/")

    def test_user_add(self):
        """Visit add page as logged in user"""
        with self.assertTemplateUsed("cyano/edit.html"):
            self.assertOK("Species/add/")

        self.assertForbidden(SPECIES + "/" + MODEL + "/add/")
    
    def test_user_add_permission(self):
        """Visit add page as logged in user with permission"""
        self.doLogin().add_allow_permission(self.getSpecies(), "WRITE_NORMAL")

        with self.assertTemplateUsed("cyano/edit.html"):
            self.assertOK("Species/add/")
            self.assertOK(SPECIES + "/" + MODEL + "/add/")
    
    def _test_editdelete_user(self, s):
        self.assertNotFound("Species/%s/" % (s))

        self.assertForbidden(SPECIES + "/%s/" % (s))
        self.assertForbidden(SPECIES + "/" + MODEL + "/" + ITEM + "/%s/" % (s))
    
    def test_user_edit(self):
        """Visit edit page as logged in user"""
        self._test_editdelete_user("edit")
    
    def test_user_edit_permission(self):
        """Visit edit page as logged in user with permission"""
        self.doLogin().add_allow_permission(self.getSpecies(), "WRITE_NORMAL")

        with self.assertTemplateUsed("cyano/edit.html"):
            self.assertOK(SPECIES + "/edit/")
            self.assertOK(SPECIES + "/" + MODEL + "/" + ITEM + "/edit/")
    
    def test_user_delete(self):
        """Visit delete page as logged in user"""
        self._test_editdelete_user("delete")
    
    def test_user_delete_permission(self):
        """Visit delete page as logged in user with permission"""
        self.doLogin().add_allow_permission(self.getSpecies(), "WRITE_DELETE")

        with self.assertTemplateUsed("cyano/delete.html"):
            self.assertOK(SPECIES + "/delete/")
            self.assertOK(SPECIES + "/" + MODEL + "/" + ITEM + "/delete/")

    def test_user_exportData(self):
        """Visit export page as logged in user"""
        self.assertForbidden(SPECIES + "/export/")
        
        self.doLogin().add_allow_permission(self.getSpecies(), "READ_NORMAL")
        
        with self.assertTemplateUsed("cyano/exportDataForm.html"):
            self.assertOK(SPECIES + "/export/")

    def test_user_import(self):
        """Visit import page as logged in user"""
        with self.assertTemplateUsed("cyano/importDataForm.html"):
            self.assertOK("/import/data/")
            self.assertOK(SPECIES + "/import/data/")
        
        with self.assertTemplateUsed("cyano/importSpeciesForm.html"):
            self.assertOK("/import/species/")
            self.assertOK(SPECIES + "/import/species/")

    def test_user_history(self):
        """Visit history pages as logged in user"""
        self.assertForbidden(SPECIES + "/history/")
        self.assertForbidden(SPECIES + "/" + MODEL + "/history/")
        self.assertForbidden(SPECIES + "/" + MODEL + "/" + ITEM + "/history/")
        
        self.doLogin().add_allow_permission(self.getSpecies(), "READ_HISTORY")
        
        with self.assertTemplateUsed("cyano/history.html"):
            self.assertOK(SPECIES + "/history/")
            self.assertOK(SPECIES + "/" + MODEL + "/history/")
            self.assertOK(SPECIES + "/" + MODEL + "/" + ITEM + "/history/")
    
    def test_user_history_detail(self):
        """Visit history detail page as logged in user"""
        self.assertForbidden(SPECIES + "/" + MODEL + "/" + ITEM + "/history/1/")
        
        self.doLogin().add_allow_permission(self.getSpecies(), "READ_HISTORY")
        
        with self.assertTemplateUsed("cyano/history_detail.html"):
            self.assertOK(SPECIES + "/" + MODEL + "/" + ITEM + "/history/1/")

    def test_user_permission(self):
        """Visit permission view page as logged in user"""    
        self.assertForbidden(SPECIES + "/permission/")
        self.assertForbidden(SPECIES + "/" + MODEL + "/" + ITEM + "/permission/")
        
        self.doLogin().add_allow_permission(self.getSpecies(), "READ_PERMISSION")
        
        with self.assertTemplateUsed("cyano/permission.html"):
            self.assertOK(SPECIES + "/permission/")
            self.assertOK(SPECIES + "/" + MODEL + "/" + ITEM + "/permission/")

    def test_user_permission_edit(self):
        """Visit permission edit page as logged in user"""
        # Special case: Matches species_wid/edit
        self.assertNotFound("permission/edit/")
        # Special case: Matches species_wid/model/edit
        self.assertNotFound(SPECIES + "/" + MODEL + "/permission/edit/")
        
        self.assertForbidden(SPECIES + "/permission/edit/")
        self.assertForbidden(SPECIES + "/" + MODEL + "/" + ITEM + "/permission/edit/")
        
        self.doLogin().add_allow_permission(self.getSpecies(), "WRITE_PERMISSION")
        
        self.assertForbidden(SPECIES + "/permission/edit/")
        self.assertForbidden(SPECIES + "/" + MODEL + "/" + ITEM + "/permission/edit/")
        
        self.doLogin().add_allow_permission(self.getSpecies(), "READ_PERMISSION")
        
        with self.assertTemplateUsed("cyano/permission_edit.html"):
            self.assertOK(SPECIES + "/permission/edit/")
            self.assertOK(SPECIES + "/" + MODEL + "/" + ITEM + "/permission/edit/")

    def test_user_sitemap(self):
        """Visit sitemaps as logged in user"""
        self._test_sitemap()

class CyanoTest(CyanoBaseTest):
    def test_worker_online(self):
        """Tests if a celery worker is online"""
        
        from celery.task.control import inspect
        
        insp = inspect()
        
        try:
            res = insp.reserved()
        except Exception:
            print "RabbitMQ offline?"
            raise

        self.assertIsNotNone(res, "No worker online")
    
    def test_is_admin(self):
        self.doLogin()
        self.assertForbidden(SPECIES + "/")
        
        self.assertTrue(self.doAdminLogin().is_admin())

        with self.assertTemplateUsed("cyano/species.html"): 
            self.assertOK(SPECIES + "/")
    
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

