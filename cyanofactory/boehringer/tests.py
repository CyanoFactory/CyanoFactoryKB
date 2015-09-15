"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""
import json

from cyano.tests import CyanoBaseTest
from django.core import management
import boehringer.models as models

ACCESS_PERM = "access_boehringer"

class BoehringerStandaloneGuestTest(CyanoBaseTest):
    """Tests standalone feature under /boehringer/ as guest"""
    def test_boehringer_data_loaded(self):
        msg = "Data not imported\nRun \"manage.py loaddata boehringer/fixtures/boehringer.json\""

        self.assertEqual(models.BioMolecule.objects.count(), 2565, msg=msg)
        self.assertEqual(models.Enzyme.objects.count(), 1016, msg=msg)
        self.assertEqual(models.Metabolite.objects.count(), 1549, msg=msg)
        self.assertEqual(models.Color.objects.count(), 3, msg=msg)
    
    def test_main_page_guest(self):
        """Main page as guest"""
        self.assertPermissionRequired("boehringer/", ACCESS_PERM)

    def test_legacy_page_guest(self):
        """Main page as guest"""
        self.assertPermissionRequired("boehringer/legacy/", ACCESS_PERM)

    def test_main_page_guest_permission(self):
        self.doGuestLogin()
        self.user.profile.assign_perm(ACCESS_PERM)

        with self.assertTemplateUsed("boehringer/index.html"):
            self.assertOK("boehringer/")

    def test_legacy_page_guest_permission(self):
        self.doGuestLogin()
        self.user.profile.assign_perm(ACCESS_PERM)

        with self.assertTemplateUsed("boehringer/legacy.html"):
            self.assertOK("boehringer/legacy/")

class BoehringerStandaloneTest(BoehringerStandaloneGuestTest):
    """Tests standalone feature under /boehringer/"""

    def setUp(self):
        self.doLogin()
    
    def test_main_page(self):
        """Main page loads"""
        self.user.profile.assign_perm(ACCESS_PERM)

        with self.assertTemplateUsed("boehringer/index.html"):
            self.assertOK("boehringer/")
    
    def test_legacy_page(self):
        """Legacy page loads"""
        self.user.profile.assign_perm(ACCESS_PERM)

        with self.assertTemplateUsed("boehringer/legacy.html"):
            self.assertOK("boehringer/legacy/")

class BoehringerStandaloneTestWithPermission(CyanoBaseTest):
    def setUp(self):
        self.doLogin()
        self.user.profile.assign_perm(ACCESS_PERM)

    def test_default_search(self):
        """default search results are correct"""
        with self.assertTemplateUsed("boehringer/index.html"):
            resp = self.assertOK("boehringer/")
            self.assertContains(resp, "4 hits - 1 miss")
            self.assertContains(resp, "L-Ascorbate")

    def test_ec_search(self):
        """EC search results are correct"""
        with self.assertTemplateUsed("boehringer/index.html"):
            resp = self.assertPOSTOK("boehringer/", data = {"items": "1.1.1.1"})
            self.assertContains(resp, "1 hit - 0 misses")
            self.assertContains(resp, "Alcohol dehydrogenase")
            self.assertContains(resp, "Aldehyde dehydrogenase")

    def test_metabolite_search(self):
        """metabolite search results are correct"""
        with self.assertTemplateUsed("boehringer/index.html"):
            resp = self.assertPOSTOK("boehringer/", data = {"items": "Atp XXXXX"})
            self.assertContains(resp, "1 hit - 1 miss")
            self.assertContains(resp, "Phosphoribosyl-ATP")
            resp = self.assertPOSTOK("boehringer/", data = {"items": "ATP XXXXX"})
            self.assertContains(resp, "1 hit - 1 miss")
            self.assertContains(resp, "Phosphoribosyl-ATP")
            resp = self.assertPOSTOK("boehringer/", data = {"items": "Atp#green XXXXX#blue"})
            self.assertContains(resp, "1 hit - 1 miss")
            self.assertContains(resp, "Phosphoribosyl-ATP")

    def test_ec_metabolite_search(self):
        """EC and metabolite search results are correct"""
        with self.assertTemplateUsed("boehringer/index.html"):
            resp = self.assertPOSTOK("boehringer/", data = {"items": "Atp XXXXX 1.1.1.1 9.9.9.9"})
            self.assertContains(resp, "1 hit - 1 miss", count = 2)
            self.assertContains(resp, "Phosphoribosyl-ATP")
            self.assertContains(resp, "Alcohol dehydrogenase")
            self.assertContains(resp, "Aldehyde dehydrogenase")

    def test_ec_metabolite_whitespace_search(self):
        """EC and metabolite search results with different whitespace are correct"""
        with self.assertTemplateUsed("boehringer/index.html"):
            resp = self.assertPOSTOK("boehringer/", data = {"items": "Atp\rXXXXX\n\n1.1.1.1\t9.9.9.9"})
            self.assertContains(resp, "1 hit - 1 miss", count = 2)
            self.assertContains(resp, "Phosphoribosyl-ATP")
            self.assertContains(resp, "Alcohol dehydrogenase")
            self.assertContains(resp, "Aldehyde dehydrogenase")

class BoehringerUsabilityTestUser(CyanoBaseTest):
    def setUp(self):
        self.doLogin()
        self.user.profile.assign_perm(ACCESS_PERM)

    def test_ajax_ok(self):
        self.assertPOSTBadRequest("boehringer/ajax/", ajax=True)

        self.assertPOSTBadRequest("boehringer/ajax/?op=invalid", ajax=True)

    def test_ajax_save(self):
        self.assertPOSTBadRequest("boehringer/ajax/?op=load&pk=1", ajax=True)

        self.assertPOSTBadRequest("boehringer/ajax/", data={"op":"save", "name": "aaa", "query": ""}, ajax=True)
        self.assertPOSTBadRequest("boehringer/ajax/", data={"op":"save", "name": "", "query": "aaa"}, ajax=True)

        res = self.assertPOSTOK("boehringer/ajax/", data={"op": "save", "name": "test", "query": "aaa"}, ajax=True)
        j = json.loads(res.content.decode())
        self.assertEquals(j["name"], "test")
        self.assertEquals(j["created"], True)

        res = self.assertPOSTOK("boehringer/ajax/", data={"op": "save", "name": "test", "query": "aaa"}, ajax=True)
        j = json.loads(res.content.decode())
        self.assertEquals(j["created"], False)

    def test_ajax_load(self):
        res = self.assertPOSTOK("boehringer/ajax/", data={"op": "save", "name": "test", "query": "aaa"}, ajax=True)
        j = json.loads(res.content.decode())

        res = self.assertPOSTOK("boehringer/ajax/", data={"op": "load", "pk": j["pk"]}, ajax=True)
        j = json.loads(res.content.decode())

        self.assertEquals(j["name"], "test")
        self.assertEquals(j["query"], "aaa")

    def test_ajax_delete(self):
        res = self.assertPOSTOK("boehringer/ajax/", data={"op": "save", "name": "test", "query": "aaa"}, ajax=True)
        j = json.loads(res.content.decode())

        self.assertPOSTOK("boehringer/ajax/", data={"op": "delete", "pk": j["pk"]}, ajax=True)

        # Deleted now
        res = self.assertPOSTBadRequest("boehringer/ajax/", data={"op": "load", "pk": j["pk"]}, ajax=True)