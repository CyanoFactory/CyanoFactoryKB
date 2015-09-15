"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""
import json

from cyano.tests import CyanoBaseTest
import kegg.models as models

ACCESS_PERM = "access_kegg"


class KeggBaseTest(CyanoBaseTest):
    """Tests standalone feature under /boehringer/ as guest"""
    def test_kegg_data_loaded(self):
        msg = "Data not imported\nRun \"manage.py import_kegg_data\""

        self.assertEqual(models.EcNumber.objects.count(), 2820, msg=msg)
        self.assertEqual(models.Map.objects.count(), 206, msg=msg)

    def test_main_page_guest(self):
        """Main page as guest"""
        self.assertPermissionRequired("kegg/", ACCESS_PERM)

    def test_main_page_guest_permission(self):
        """Main page as guest with permission"""
        self.doGuestLogin()
        self.user.profile.assign_perm(ACCESS_PERM)

        with self.assertTemplateUsed("kegg/index.html"):
            self.assertOK("kegg/")


class KeggUsabilityTest(CyanoBaseTest):
    def setUp(self):
        self.doGuestLogin()
        self.user.profile.assign_perm(ACCESS_PERM)

    def test_access_map(self):
        """Access a map"""
        with self.assertTemplateUsed("kegg/map.html"):
            self.assertOK("kegg/map00010/")

        with self.assertTemplateUsed("kegg/map.html"):
            self.assertOK("kegg/map01120/")

    def text_ajax_ok(self):
        self.assertBadRequest("kegg/ajax/")

        # Guests can't access ajax
        self.assertLoginRequired("kegg/ajax/", ajax=True)


class KeggUsabilityTestUser(CyanoBaseTest):
    def setUp(self):
        self.doLogin()
        self.user.profile.assign_perm(ACCESS_PERM)

    def test_access_map(self):
        """Access a map"""
        with self.assertTemplateUsed("kegg/map.html"):
            self.assertOK("kegg/map00010/")

        with self.assertTemplateUsed("kegg/map.html"):
            self.assertOK("kegg/map01120/")

    def test_ajax_ok(self):
        self.assertPOSTBadRequest("kegg/ajax/", ajax=True)

        self.assertPOSTBadRequest("kegg/ajax/?op=invalid", ajax=True)

    def test_ajax_save(self):
        self.assertPOSTBadRequest("kegg/ajax/?op=load&pk=1", ajax=True)

        self.assertPOSTBadRequest("kegg/ajax/", data={"op":"save", "name": "aaa", "query": ""}, ajax=True)
        self.assertPOSTBadRequest("kegg/ajax/", data={"op":"save", "name": "", "query": "aaa"}, ajax=True)

        res = self.assertPOSTOK("kegg/ajax/", data={"op": "save", "name": "test", "query": "aaa"}, ajax=True)
        j = json.loads(res.content.decode())
        self.assertEquals(j["name"], "test")
        self.assertEquals(j["created"], True)

        res = self.assertPOSTOK("kegg/ajax/", data={"op": "save", "name": "test", "query": "aaa"}, ajax=True)
        j = json.loads(res.content.decode())
        self.assertEquals(j["created"], False)

    def test_ajax_load(self):
        res = self.assertPOSTOK("kegg/ajax/", data={"op": "save", "name": "test", "query": "aaa"}, ajax=True)
        j = json.loads(res.content.decode())

        res = self.assertPOSTOK("kegg/ajax/", data={"op": "load", "pk": j["pk"]}, ajax=True)
        j = json.loads(res.content.decode())

        self.assertEquals(j["name"], "test")
        self.assertEquals(j["query"], "aaa")

    def test_ajax_delete(self):
        res = self.assertPOSTOK("kegg/ajax/", data={"op": "save", "name": "test", "query": "aaa"}, ajax=True)
        j = json.loads(res.content.decode())

        self.assertPOSTOK("kegg/ajax/", data={"op": "delete", "pk": j["pk"]}, ajax=True)

        # Deleted now
        res = self.assertPOSTBadRequest("kegg/ajax/", data={"op": "load", "pk": j["pk"]}, ajax=True)