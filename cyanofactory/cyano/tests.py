'''
CyanoFactory knowledge base tests

Author: Jonathan Karr, jkarr@stanford.edu
Affiliation: Covert Lab, Department of Bioengineering, Stanford University
Last updated: 2012-07-17

Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
'''
from django.core.urlresolvers import reverse

from django.test import TestCase, Client
from cyano.helpers import getModelsMetadata

import cyano.models as cmodels
from cyano.models import PermissionEnum as P
from django.core.exceptions import ObjectDoesNotExist
from guardian.shortcuts import assign_perm
from django.contrib.auth.models import User, Group

import settings
from cyano.forms import ImportSpeciesForm

SPECIES = "Synechocystis"
SPECIES2 = "Ecoli"
MODEL = "Gene"
ITEM = "sll7106"

#def setup():
    ## Permissions
    #management.call_command('loaddata', 'metadata.json')
    ## Create users and groups
    #management.call_command('autocreateinitial')
    ## Create meta data tables
    #management.call_command('create_meta')
    ## Create simple entry
    #g = Genbank(wid = SPECIES,
    #        user = "management",
    #        reason = "Unit Testing",
    #        chromosome = "UnitTest-Chromosome",
    #        name = "Chromosome for Unit Test")
    #with open("../sequences/testcase.gb", "r") as f:
    #    g.parse(f)
    #g.apply()
    #
    #species = cmodels.Species(wid = "Other-Species", name = "Unit Test species 2")
    #rev = cmodels.RevisionDetail(
    #    user = cmodels.UserProfile.objects.get(user__username = "management"),
    #    reason = "Unit Testing")
    #rev.save()
    #species.save(rev)

class CyanoBaseTest(TestCase):
    def __init__(self, *args, **kwargs):
        self.user = None

        super(TestCase, self).__init__(*args, **kwargs)

    def setUp(self):
        self.client = Client()
    
    def doLogin(self):
        self.client.login(username="gabriel", password="aaa")
        self.user = User.objects.get(username="gabriel")
    
    def doLogout(self):
        self.client.logout()
    
    def doAdminLogin(self):
        self.client.login(username="admin", password="aaa")
        self.user = User.objects.get(username="admin")

    def doGuestLogin(self):
        """Fake login returning the guest handle"""
        self.user = User.objects.get(username="guest")
    
    def getSpecies(self, second = False):
        return cmodels.Species.objects.for_wid(SPECIES2 if second else SPECIES)
    
    def getItem(self):
        return getattr(cmodels, MODEL).objects.for_species(SPECIES).for_wid(ITEM)

    def checkMany(self, li, assert_func):
        for l in li:
            assert_func(l)

    # via http://djangosnippets.org/snippets/137/
    def GET(self, url, params={}, status=200, mimetype="text/html", follow=False, ajax=False):
        """Get a URL and require a specific status code before proceeding"""
        if url[0] != "/":
            url = "/" + url
        response = self.client.get(url, params, follow = follow, HTTP_X_REQUESTED_WITH='XMLHttpRequest' if ajax else "")

        if response.status_code != status:
            # print note that login was required when a test fails
            if "location" in response._headers:
                if "login" in response._headers["location"][1]:
                    print("Login required.")

        self.assertEqual(response.status_code, status, "\nURL: {}\nStatus {} (expected: {})".format(url, response.status_code, status))
            
        self.assertTrue(response["content-type"].startswith(mimetype), "Mimetype {} (expected: {})".format(response["content-type"], mimetype))     
        return response

    def POST(self, url, params = {}, status=200, mimetype="text/html", ajax=False):
        """Make a POST and require a specific status code before proceeding"""
        if url[0] != "/":
            url = "/" + url
        response = self.client.post(url, params, HTTP_X_REQUESTED_WITH='XMLHttpRequest' if ajax else "")
        self.assertEqual(response.status_code, status, "\nURL: {}\nStatus {} (expected: {})".format(url, response.status_code, status))
        self.assertTrue(response["content-type"].startswith(mimetype), "Mimetype {} (expected: {})".format(response["content-type"], mimetype))
        return response

    def assertOK(self, url, data={}, follow=False, mimetype="text/html", ajax=False):
        """Fail if result code is not 200"""
        return self.GET(url, params=data, status=200, follow=follow, mimetype=mimetype, ajax=ajax)
    
    def assertPOSTOK(self, url, data={}, follow=False, ajax=False):
        """Fail if result code is not 200"""
        mimetype = "application/json" if ajax else "text/html"

        return self.POST(url, params=data, status=200, ajax=ajax, mimetype=mimetype)

    def assertNotFound(self, url, follow = False):
        """Fail if found (code != 404)"""
        return self.GET(url, status=404, follow=follow)
    
    def assertForbidden(self, url, follow = False):
        """Fail if not forbidden (code 403)"""
        return self.GET(url, status=403, follow=follow)

    def assertBadRequest(self, url, follow = False, ajax=False, data={}):
        """Fail if not bad request (code 400)"""
        mimetype = "application/json" if ajax else "text/html"

        return self.GET(url, status=400, follow=follow, ajax=ajax, params=data, mimetype=mimetype)

    def assertPOSTBadRequest(self, url, follow=False, ajax=False, data={}):
        """Fail if not bad request (code 400)"""
        mimetype = "application/json" if ajax else "text/html"

        return self.POST(url, status=400, ajax=ajax, params=data, mimetype=mimetype)
    
    def assertRedirect(self, url, target_url=None, target_status=200):
        """Fail if a url does not redirect (code 302)"""
        response = self.GET(url, status = 302)
        if target_url:
            if target_url[0] != "/":
                target_url = "/" + target_url
            return self.assertRedirects(response, target_url, 302, target_status)
        return response
    
    def assertRedirectPermanent(self, url, target_url=None, target_status=200):
        """Fail if a url does not permanently redirect (code 301)"""
        response = self.GET(url, status = 301)
        if target_url:
            if target_url[0] != "/":
                target_url = "/" + target_url
            return self.assertRedirects(response, target_url, 301, target_status)
        return response
    
    def assertMimeTypeText(self, url, status=200, follow=False):
        """Fail if mimetype is not text"""
        return self.GET(url, status = status, mimetype = "text/plain", follow=follow)
    
    def assertMimeTypeJson(self, url, status=200, follow = False):
        """Fail if mimetype is not json"""
        return self.GET(url, status = status, mimetype = "application/json", follow = follow)
    
    def assertMimeTypeXml(self, url, status=200, follow = False):
        """Fail if mimetype is not XML"""
        return self.GET(url, status = status, mimetype = "application/xml", follow = follow)

    def assertLoginRequired(self, url, ajax=False):
        """Tests if a view redirects to login page"""
        # Not for testing the decorator, only if the decorator is set

        with self.assertTemplateUsed("registration/login.html"):
            self.assertOK(url, follow=True, ajax=ajax)

    def assertPermissionRequired(self, url, required_perm):
        """Tests if a view redirects to login page and passes a no permission status text"""
        with self.assertTemplateUsed("registration/login.html"):
            response = self.assertForbidden(url)
            self.assertEqual(response.context_data.get("required_perm", "").lower(),
                              required_perm.lower())

    def assignPerm(self, which, obj=None):
        self.user.profile.assign_perm(which, obj)

class CyanoBasicTest(CyanoBaseTest):
    def test_logic_works(self):
        self.assertFalse(True is False, "Hardware error")

    def test_append_slash(self):
        """Test if APPEND_SLASH setting is enabled"""
        self.doLogin()
        
        with self.assertTemplateUsed("registration/login.html"): 
            self.assertRedirectPermanent("login", "login/")

    def test_worker_online(self):
        """Tests if a celery worker is online"""
        
        from celery.task.control import inspect
        
        insp = inspect()
        
        try:
            res = insp.reserved()
        except Exception:
            print("RabbitMQ offline?")
            raise

        self.assertIsNotNone(res, "No worker online")
    
    def test_is_admin(self):
        self.doLogin()
        self.assertForbidden(SPECIES + "/")

        self.doAdminLogin()
        self.assertTrue(self.user.is_superuser)

        with self.assertTemplateUsed("cyano/species.html"): 
            self.assertOK(SPECIES + "/")

class CyanoTest(CyanoBaseTest):
    def test_permissions_species(self):
        self.assertPermissionRequired(SPECIES + "/", P.ACCESS_SPECIES)
        
        self.doGuestLogin()
        self.user.profile.assign_perm(P.ACCESS_SPECIES)

        self.assertPermissionRequired(SPECIES + "/", P.READ_NORMAL)

        entry = cmodels.Species.objects.for_wid(SPECIES)
        self.user.profile.assign_perm(P.READ_NORMAL, entry)

        with self.assertTemplateUsed("cyano/species.html"): 
            self.assertOK(SPECIES + "/")
    
    #def test_empty_permission_deleted(self):
    #    user = self.doLogin()
    #    entry = cmodels.Species.objects.for_wid(SPECIES)
    #
    #    with self.assertRaises(ObjectDoesNotExist):
    #        cmodels.UserPermission.objects.get(user = user, entry = entry)
    #
    #    user.add_allow_permission(entry, "READ_NORMAL")
   # #    cmodels.UserPermission.objects.get(user = user, entry = entry)
    #    user.delete_allow_permission(entry, "READ_NORMAL")
    #
    #    with self.assertRaises(ObjectDoesNotExist):
    #        cmodels.UserPermission.objects.get(user = user, entry = entry)
    
    def test_species_not_found(self):
        self.doGuestLogin()
        self.user.profile.assign_perm(P.ACCESS_SPECIES)

        with self.assertTemplateUsed("cyano/error.html"):
            self.assertNotFound("/Wrong-Species/")

class CyanoUserTestBase(CyanoBaseTest):
    def _test_index(self):
        with self.assertTemplateUsed("cyano/index.html"): 
            self.assertOK("/")

    def _test_loginout(self, s, tmpl):
        with self.assertTemplateUsed("registration/%s.html" % (tmpl)):
            self.assertOK("%s/" % (s))
            self.assertOK(SPECIES + "/%s/" % (s))

        self.assertNotFound("%s/invalid/" % (s))
        self.assertNotFound(SPECIES + "/%s/invalid/" % s)
        self.assertNotFound(SPECIES + "/model/%s/" % s)

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


class CyanoUserTestAccessSpeciesPerm(CyanoUserTestBase):
    """Tests that already set ACCESS_SPECIES permission"""

    def setUp(self):
        """Give Guest ACCESS_SPECIES permission, this check is used anywhere"""

        super(CyanoUserTestAccessSpeciesPerm, self).setUp()
        self.doGuestLogin()
        self.user.profile.assign_perm(P.ACCESS_SPECIES)

    def test_index(self):
        """Visits the main page"""
        self._test_index()

    def test_login(self):
        """Visit login page"""
        self._test_loginout("login", "login")
    
    def test_logout(self):
        """Visit logout page"""
        self._test_loginout("logout", "logged_out")

    def _test_basic_uri_stuff(self, s):
        self.assertLoginRequired("%s/" % s)
        self.assertLoginRequired(SPECIES + "/%s/" % s)
        
        self.assertNotFound("%s/invalid/" % s)
        self.assertNotFound(SPECIES + "/%s/invalid/" % s)
        self.assertNotFound(SPECIES + "/%s/user/" % s)

    def test_users(self):
        """Visit users page"""
        self._test_basic_uri_stuff("user")
    
    def test_user(self):
        """Visit user profile pages"""
        self._test_basic_uri_stuff("user/admin")
        self._test_basic_uri_stuff("user/invalid")
    
    def test_job(self):
        """Visit jobs page"""
        self._test_basic_uri_stuff("jobs")
    
    def test_search_google_on(self):
        """Tests search functions (google search enabled)"""
        settings.GOOGLE_SEARCH_ENABLED = True

        self._test_search(google=True)

    def test_search_google_off(self):
        """Tests search functions (google search disabled)"""
        settings.GOOGLE_SEARCH_ENABLED = False

        self._test_search(google=False)

    def test_species(self):
        """Visit species page"""
        self.assertPermissionRequired(SPECIES + "/", P.READ_NORMAL)
        self.assignPerm(P.READ_NORMAL, self.getSpecies())

        with self.assertTemplateUsed("cyano/species.html"):
            self.assertOK(SPECIES + "/")
    
    def test_listing(self):
        """Visit listing page"""
        url = SPECIES + "/" + MODEL + "/"
        
        self.assertPermissionRequired(url, P.READ_NORMAL)
        self.assertNotFound(SPECIES + "/" + "Invalid" + "/")

        self.assignPerm(P.READ_NORMAL, self.getSpecies())
        
        self.assertNotFound(SPECIES + "/" + "Invalid" + "/")
        
        with self.assertTemplateUsed("cyano/list.html"):
            self.assertOK(url)
        
        self.assertNotFound(SPECIES + "/Species/")
        self.assertNotFound(SPECIES + "/Entry/")

    def test_listing_many(self):
        """Visit different listing pages"""
        self.assignPerm(P.READ_NORMAL, self.getSpecies())

        for key in getModelsMetadata(cmodels.SpeciesComponent).keys():
            with self.assertTemplateUsed("cyano/list.html"):
                self.assertOK(reverse("cyano.views.listing", kwargs=dict(species_wid=SPECIES, model_type=key)))

    def test_detail(self):
        """Visit detail page"""
        url = SPECIES + "/" + MODEL + "/" + ITEM + "/"
        url2 = SPECIES2 + "/" + MODEL + "/" + ITEM + "/"
        self.assertPermissionRequired(url, P.READ_NORMAL)
        self.assertNotFound(url2)

        self.assertNotFound(SPECIES + "/" + "Invalid" + "/" + ITEM + "/")
        self.assertNotFound(SPECIES + "/" + MODEL + "/" + "Invalid" + "/")

        self.assignPerm(P.READ_NORMAL, self.getSpecies())

        self.assertNotFound(SPECIES + "/" + "Invalid" + "/" + ITEM + "/")
        self.assertNotFound(SPECIES + "/" + MODEL + "/" + "Invalid" + "/")
        
        with self.assertTemplateUsed("cyano/detail.html"):
            self.assertOK(url)

    def test_add(self):
        """Visit add page"""

        # Special case for species
        self.assertNotFound("Species/")
        self.assertPermissionRequired("Species/add/", P.ACCESS_SPECIES)

        self.assertPermissionRequired(SPECIES + "/" + MODEL + "/add/", P.ACCESS_SPECIES)
        self.assertPermissionRequired(SPECIES + "/add/", P.ACCESS_SPECIES)
        self.assertPermissionRequired(SPECIES + "/add/invalid/")

        self.assertNotFound(SPECIES + "/" + MODEL + "/" + ITEM + "/add/")
        self.assertNotFound(SPECIES + "/" + MODEL + "/" + ITEM + "/add/invalid/")
        self.assertNotFound(SPECIES + "/" + MODEL + "/add/invalid/")
    
    def test_add_permission(self):
        """Visit add page with permission"""
        self.doGuestLogin()
        self.user.profile.assign_perm(P.ACCESS_SPECIES)

        self.assertPermissionRequired(SPECIES + "/" + MODEL + "/add/", P.WRITE_NORMAL)
        self.user.profile.assign_perm(P.WRITE_NORMAL, self.getSpecies())

        with self.assertTemplateUsed("cyano/edit.html"):
            self.assertOK(SPECIES + "/" + MODEL + "/add/")
    
    def _test_editdelete(self, s, p):
        self.assertNotFound("%s/" % s)
        self.assertNotFound("Species/%s/")

        self.assertNotFound(SPECIES + "/" + MODEL + "/%s/" % s)
        self.assertNotFound(SPECIES + "/" + MODEL + "/" + ITEM + "/%s/invalid/" % s)
        self.assertNotFound(SPECIES + "/%s/invalid/" % s)
        self.assertNotFound(SPECIES + "/" + MODEL + "/%s/invalid/" % s)
        
        self.assertPermissionRequired(SPECIES + "/%s/" % s, p)
        self.assertPermissionRequired(SPECIES + "/" + MODEL + "/" + ITEM + "/%s/" % s, p)
    
    def test_edit(self):
        """Visit edit page"""
        self._test_editdelete("edit", P.WRITE_NORMAL)
    
    def test_edit_permission(self):
        """Visit edit page with permission"""
        self.assignPerm(P.WRITE_NORMAL, self.getSpecies())

        with self.assertTemplateUsed("cyano/edit.html"):
            self.assertOK(SPECIES + "/edit/")
    
    def test_delete(self):
        """Visit delete page"""
        self._test_editdelete("delete", P.WRITE_DELETE)
    
    def test_delete_permission(self):
        """Visit delete page with permission"""
        self.assignPerm(P.WRITE_DELETE, self.getSpecies())

        with self.assertTemplateUsed("cyano/delete.html"):
            self.assertOK(SPECIES + "/delete/")

    def test_exportData(self):
        """Visit export page"""
        self.assertPermissionRequired(SPECIES + "/export/", P.READ_NORMAL)
        
        self.assignPerm(P.READ_NORMAL, self.getSpecies())
        
        with self.assertTemplateUsed("cyano/exportDataForm.html"):
            self.assertOK(SPECIES + "/export/")
        
        self.assertNotFound("export/")
        self.assertNotFound("invalid/export/")
        self.assertNotFound(SPECIES + "/model/export/")

    def test_import(self):
        """Visit import page"""
        self.assertNotFound("import/")
        
        self.assertNotFound("/import/data/")
        self.assertNotFound("/import/species/")
        
        self.assertPermissionRequired(SPECIES + "/import/data/", P.WRITE_NORMAL)
        self.assertPermissionRequired(SPECIES + "/import/species/", P.WRITE_NORMAL)
        
        self.assertNotFound("/import/data/invalid/")
        self.assertNotFound("/import/species/invalid/")
        self.assertNotFound(SPECIES + "/import/data/invalid/")
        self.assertNotFound(SPECIES + "/import/species/invalid/")

    def test_history(self):
        """Visit history pages"""
        self.assertNotFound("history/")
        
        self.assertPermissionRequired(SPECIES + "/history/", P.READ_HISTORY)
        self.assertPermissionRequired(SPECIES + "/" + MODEL + "/history/", P.READ_HISTORY)
        self.assertPermissionRequired(SPECIES + "/" + MODEL + "/" + ITEM + "/history/", P.READ_HISTORY)
        
        self.assignPerm(P.READ_HISTORY, self.getSpecies())
        
        with self.assertTemplateUsed("cyano/history.html"):
            self.assertOK(SPECIES + "/history/")
            self.assertOK(SPECIES + "/" + MODEL + "/history/")
            self.assertOK(SPECIES + "/" + MODEL + "/" + ITEM + "/history/")
    
    def test_history_detail(self):
        """Visit history detail page"""
        self.assertNotFound(SPECIES + "/" + MODEL + "/" + ITEM + "/history/a/")
        self.assertNotFound(SPECIES + "/" + MODEL + "/" + ITEM + "/history/1/a/")
        
        self.assertPermissionRequired(SPECIES + "/" + MODEL + "/" + ITEM + "/history/1/", P.READ_HISTORY)
        
        self.assignPerm(P.READ_HISTORY, self.getSpecies())
        
        with self.assertTemplateUsed("cyano/history_detail.html"):
            self.assertOK(SPECIES + "/" + MODEL + "/" + ITEM + "/history/3/")

    def test_permission(self):
        """Visit permission view page"""
        #self.assertRedirect("permission/", "login/?next={}".format(reverse('cyano.views.global_permission')))
        
        self.assertPermissionRequired(SPECIES + "/permission/", P.READ_PERMISSION)
        self.assertNotFound(SPECIES + "/" + MODEL + "/permission/")
        self.assertPermissionRequired(SPECIES + "/" + MODEL + "/" + ITEM + "/permission/", P.READ_PERMISSION)
        
        self.assignPerm(P.READ_PERMISSION, self.getSpecies())
        
        with self.assertTemplateUsed("cyano/global_permission.html"):
            self.assertOK(SPECIES + "/permission/")
            self.assertOK(SPECIES + "/" + MODEL + "/" + ITEM + "/permission/")

    def test_permission_edit(self):
        """Visit permission edit page"""
        # Global permissions have no edit url
        self.assertNotFound("permission/edit/")
        # model has no permissions
        self.assertNotFound(SPECIES + "/" + MODEL + "/permission/edit/")
        
        self.assertPermissionRequired(SPECIES + "/permission/edit/", P.WRITE_PERMISSION)
        self.assertPermissionRequired(SPECIES + "/" + MODEL + "/" + ITEM + "/permission/edit/", P.WRITE_PERMISSION)

        self.assignPerm(P.WRITE_PERMISSION, self.getSpecies())

        self.assertOK(SPECIES + "/permission/edit/")
        self.assertOK(SPECIES + "/" + MODEL + "/" + ITEM + "/permission/edit/")

    def test_sitemap(self):
        """Visit sitemaps"""
        self.assertNotFound("sitemap.xml/")
        self.assertNotFound("sitemap_toplevel.xml/")
        self.assertNotFound("sitemap/%s.xml/" % (SPECIES))

        with self.assertTemplateUsed("cyano/sitemap.xml"):
            self.assertMimeTypeXml("sitemap.xml", mimetype="application/xml")

        with self.assertTemplateUsed("cyano/sitemap_toplevel.xml"):
            self.assertMimeTypeXml("sitemap_toplevel.xml", mimetype="application/xml")

        self.assertPermissionRequired("sitemap/%s.xml" % SPECIES, P.READ_NORMAL)
        self.assignPerm(P.READ_NORMAL, self.getSpecies())
        with self.assertTemplateUsed("cyano/sitemap_species.xml"):
            self.assertMimeTypeXml("sitemap/%s.xml" % SPECIES, mimetype="application/xml")

        with self.assertTemplateUsed("cyano/error.html"):
            self.assertNotFound("sitemap/Wrong-Species.xml")
    
    def test_license(self):
        """Visit license page"""
        with self.assertTemplateUsed("cyano/license.html"):
            self.assertOK("license/")
            self.assertOK(SPECIES + "/license/")

        self.assertNotFound("license/invalid/")
        self.assertNotFound(SPECIES + "/license/invalid/")
        self.assertNotFound(SPECIES + "/model/license/")

class CyanoUserTestReadNormalPerm(CyanoUserTestBase):
    """Tests that already set ACCESS_SPECIES and READ_NORMAL permission"""

    def setUp(self):
        """Give Guest ACCESS_SPECIES and READ_NORMAL permission, this check is used anywhere"""

        super(CyanoUserTestReadNormalPerm, self).setUp()
        self.doGuestLogin()
        self.user.profile.assign_perm(P.ACCESS_SPECIES)
        self.assignPerm(P.READ_NORMAL, self.getSpecies())

    def _test_detail_formats(self, url, fasta_fail=False, genbank_fail=False):
        with self.assertTemplateUsed("cyano/detail.html"):
            self.assertOK(url)

        # Not testing json, is rest-framework code like xml
        self.assertMimeTypeXml("/api/" + url + "?format=xml")

        def ok_func(u):
            return self.assertOK(u, mimetype="application/octet-stream")

        fasta_func = self.assertBadRequest if fasta_fail else ok_func
        genbank_func = self.assertBadRequest if genbank_fail else ok_func

        fasta_func(url + "?format=fasta")
        genbank_func(url + "?format=genbank")

    def test_detail_plasmid(self):
        self.assertNotFound(SPECIES + "/Chromosome/pSYSA/")

        self._test_detail_formats(SPECIES + "/Plasmid/pSYSA/")

    def test_detail_compartment(self):
        url = SPECIES + "/Compartment/Cytosol/"

        args = [True, True]

        self._test_detail_formats(SPECIES + "/Compartment/Cytosol/", *args)
        self._test_detail_formats(SPECIES + "/Compartment/Extra_cellular/", *args)

    def test_detail_genes(self):
        args = [False, True]

        # mRNA
        self._test_detail_formats(SPECIES + "/Gene/sll7001/", *args)

        # rRNA
        self._test_detail_formats(SPECIES + "/Gene/ST6803r01/", *args)

        # tRNA, reversed
        self._test_detail_formats(SPECIES + "/Gene/ST6803t18/", *args)

    def test_detail_metabolites(self):
        args = [True, True]

        self._test_detail_formats(SPECIES + "/Metabolite/ATP/", *args)
        self._test_detail_formats(SPECIES + "/Metabolite/ethanol/", *args)

    def test_detail_pathways(self):
        args = [True, True]

        self._test_detail_formats(SPECIES + "/Pathway/Boehringer/", *args)
        self._test_detail_formats(SPECIES + "/Pathway/map00010/", *args)

    def test_detail_proteinmonomers(self):
        self._test_detail_formats(SPECIES + "/ProteinMonomer/sll7001_Monomer/")
        self._test_detail_formats(SPECIES + "/ProteinMonomer/sll0575_Monomer/")

    def test_detail_reactions(self):
        args = [True, True]

        self._test_detail_formats(SPECIES + "/Reaction/_1_1_1_100a_r/", *args)
        self._test_detail_formats(SPECIES + "/Reaction/Biomass_r/", *args)

    def test_detail_transcription_units(self):
        args = [False, True]

        self._test_detail_formats(SPECIES + "/TranscriptionUnit/TU_sll0006-sll0007-sll0008/", *args)
        self._test_detail_formats(SPECIES + "/TranscriptionUnit/TU_sll0661/", *args)

    def test_detail_publications(self):
        args = [True, True]

        self._test_detail_formats(SPECIES + "/PublicationReference/PUB_14686584/", *args)
        self._test_detail_formats(SPECIES + "/PublicationReference/REF_0001/", *args)


class CyanoLoginTest(CyanoUserTestBase):
    def setUp(self):
        super(CyanoLoginTest, self).setUp()
        self.doLogin()
        self.user.profile.assign_perm(P.ACCESS_SPECIES)
    
    def test_user_index(self):
        """Visits the main page as logged in user"""
        self._test_index()

    def test_user_login(self):
        """Visit login page as logged in user"""
        self._test_loginout("login", "login")
    
    def test_user_logout(self):
        """Visit logout page as logged in user"""
        self._test_loginout("logout", "logged_out")

    def test_user_users(self):
        """Visit users page as logged in user"""
        with self.assertTemplateUsed("cyano/users.html"):
            self.assertOK("user/")
            self.assertOK(SPECIES + "/user/")
    
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

class CyanoCreateTest(CyanoBaseTest):
    def setUp(self):
        super(CyanoCreateTest, self).setUp()
        self.doGuestLogin()
        self.user.profile.assign_perm(P.ACCESS_SPECIES)
        self.doLogin()
        self.user.profile.assign_perm(P.ACCESS_SPECIES)
    
    def test_new_species_missing(self):
        """Test new species is missing"""
        with self.assertRaises(ObjectDoesNotExist):
            cmodels.Species.objects.for_wid("new-species")

    def test_create_species(self):
        """Create a new species via the webinterface"""
        # Form validation logic tested below
        with self.assertTemplateUsed("cyano/importSpeciesForm.html"):
            data = {
                'new_species': "New Species",
                'new_wid': "new-species",
                'reason': 'unit test'
            }
            self.assertPOSTOK("/import/species/", data=data)

        cmodels.Species.objects.for_wid("new-species")
        
        self.assertNotFound("new-species/" + MODEL + "/" + ITEM + "/")
        
        # Permission test, POST created permissions
        with self.assertTemplateUsed("cyano/species.html"):
            self.assertOK("new-species/")
        
        self.doLogout()
        self.assertPermissionRequired("new-species/", P.READ_NORMAL)

    # Missing fields are not tested because all fields have required = True -> django feature
    
    def test_create_species_wid_used(self):
        """SpeciesForm Wid in use"""
        data = {
            'new_species': "New Species",
            'new_wid': SPECIES,
            'reason': 'unit test'
        }
        form = ImportSpeciesForm(data=data)
        self.assertFalse(form.is_valid())
        self.assertFalse("new_wid" in form.cleaned_data)

    def test_create_species_wid_invalid(self):
        """SpeciesForm Wid invalid"""
        data = {
            'new_species': "New Species",
            'new_wid': "bad wid",
            'reason': 'unit test'
        }
        form = ImportSpeciesForm(data=data)
        self.assertFalse(form.is_valid())
        self.assertFalse("new_wid" in form.cleaned_data)

    def test_create_species_ok(self):
        """SpeciesForm all fine"""
        data = {
            'new_species': "New Species",
            'new_wid': "good-wid",
            'reason': 'unit test'
        }
        form = ImportSpeciesForm(data=data)
        self.assertTrue(form.is_valid())
