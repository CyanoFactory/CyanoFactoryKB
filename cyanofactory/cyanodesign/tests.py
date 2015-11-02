import json
from django.core.exceptions import ObjectDoesNotExist
from PyNetMet2.metabolism import Metabolism
from cyano.tests import CyanoBaseTest
from cyanodesign.json_model import JsonModel
from cyanodesign.models import DesignTemplate, DesignModel

ACCESS_PERM = "access_cyanodesign"
MODEL_FILE  = "../sample_data/Toy_model_VLCsubject.txt"
MODEL_FILE2 = "../sample_data/iSyn811_v3-may15_v2.txt"
MODEL_FILE3 = "../sample_data/iSyn811_v3-may15_v2.xml"

class CyanoDesignGuestTest(CyanoBaseTest):
    def test_main_page_guest(self):
        """Main page as guest"""
        self.assertPermissionRequired("design/", ACCESS_PERM)

class CyanoDesignTest(CyanoDesignGuestTest):
    """Tests features under design/ with login"""

    def setUp(self):
        self.doLogin()

    def test_main_page_perm(self):
        self.user.profile.assign_perm(ACCESS_PERM)

        with self.assertTemplateUsed("cyanodesign/list.html"):
            self.assertOK("design/")

class CyanoDesignPermissionTest(CyanoBaseTest):
    """Test features under design/ with login and permission"""

    def setUp(self):
        self.doLogin()
        self.user.profile.assign_perm(ACCESS_PERM)

    def test_template_usage(self):
        """Test usage of templates, though no template was imported"""
        form_data = dict(name="Test", choice=1)
        res = self.assertPOSTOK("design/upload/2/", data=form_data, ajax=True)
        self.assertFalse(json.loads(res.content.decode("utf-8"))["success"])

    def test_model_upload(self):
        with open(MODEL_FILE, "r") as f:
            form_data = dict(name="Test", file=f)
            res = self.assertPOSTOK("design/upload/1/", data=form_data, ajax=True)
            self.assertTrue(json.loads(res.content.decode("utf-8"))["success"])

class CyanoDesignModelTest(CyanoBaseTest):
    """Test features under design/ with login, permission and reploaded templates"""

    def setUp(self):
        self.doLogin()
        self.user.profile.assign_perm(ACCESS_PERM)

        # Add templates
        m = Metabolism(MODEL_FILE)
        j = JsonModel.from_model(m)
        dt = DesignTemplate(name="Test")
        dt.description = "Test"
        dt.filename = "ToyModel.txt"
        dt.content = j.to_json()
        dt.save()
        self.pk = dt.pk

        m = Metabolism(MODEL_FILE2)
        j = JsonModel.from_model(m)
        dt = DesignTemplate(name="Test")
        dt.description = "Test"
        dt.filename = "iSyn.txt"
        dt.content = j.to_json()
        dt.save()

    def test_template_usage(self):
        form_data = dict(name="Test", choice=self.pk)
        res = self.assertPOSTOK("design/upload/2/", data=form_data, ajax=True)
        self.assertTrue(json.loads(res.content.decode("utf-8"))["success"])

        form_data = dict(name="Test", choice=self.pk+1)
        res = self.assertPOSTOK("design/upload/2/", data=form_data, ajax=True)
        self.assertTrue(json.loads(res.content.decode("utf-8"))["success"])

        form_data = dict(name="Test", choice=self.pk+2)
        res = self.assertPOSTOK("design/upload/2/", data=form_data, ajax=True)
        self.assertFalse(json.loads(res.content.decode("utf-8"))["success"])

class CyanoDesignModelTestSimulate(CyanoDesignModelTest):
    def setUp(self):
        super().setUp()

        form_data = dict(name="Test", choice=self.pk)
        res = self.assertPOSTOK("design/upload/2/", data=form_data, ajax=True)

        form_data = dict(name="Test", choice=self.pk+1)
        res = self.assertPOSTOK("design/upload/2/", data=form_data, ajax=True)

    def test_export(self):
        dm = DesignModel.objects.all()

        for m in dm:
            url = "design/export/{}/".format(m.pk)
            self.assertOK(url, mimetype="application/x-bioopt")
            self.assertOK(url, data={"format": "bioopt"}, mimetype="application/x-bioopt")
            self.assertOK(url, data={"format": "sbml"}, mimetype="application/sbml+xml")
            self.assertBadRequest(url, data={"format": "blub"})

    def test_simulation_export(self):
        import json

        dm = DesignModel.objects.all()
        data = dict(changes="[]", display="[]", auto_flux="true")

        for m in dm:
            url = "design/simulate/{}/".format(m.pk)
            objectives = json.dumps(m.get_latest_revision().content["objectives"])
            data.update(objectives=objectives)

            data.update(dict(format="json"))
            self.assertPOSTOK(url, data=data, mimetype="application/json")

            data.update(dict(format="csv"))
            self.assertPOSTOK(url, data=data, mimetype="text/csv")

            #data.update(dict(format="png"))
            #self.assertPOSTOK(url, data=data, mimetype="image/png")

            data.update(dict(format="svg"))
            self.assertPOSTOK(url, data=data, mimetype="image/svg+xml")

            data.update(dict(format="fail"))
            self.assertBadRequest(url, data=data)


class CyanoDesignFBATest(CyanoBaseTest):
    """Tests FBA in PyNetMet and the JSON Model"""
    def test_toy_model_fba(self):
        m = Metabolism(MODEL_FILE)
        flux = 3.0
        self.assertEquals(m.fba().Z, flux)
        m = JsonModel.from_model(m).to_model()
        self.assertEquals(m.fba().Z, flux)

        m.get_reaction("reac1").constraint = (0, 5)
        flux = 5.5
        self.assertEquals(m.fba().Z, flux)
        m = JsonModel.from_model(m).to_model()
        self.assertEquals(m.fba().Z, flux)

        m.get_reaction("reac1").stoic = [[1.5],[2.5]]
        flux = 9.25
        self.assertEquals(m.fba().Z, flux)
        m = JsonModel.from_model(m).to_model()
        self.assertEquals(m.fba().Z, flux)

    def test_isyn_model_bioopt(self):
        m = Metabolism(MODEL_FILE2)
        flux = 0.08940293757452249
        self.assertEquals(m.fba().Z, flux)
        m = JsonModel.from_model(m).to_model()
        self.assertEquals(m.fba().Z, flux)

    def test_isyn_model_sbml(self):
        m = Metabolism(MODEL_FILE3)
        flux = 0.08940293757452249
        self.assertEquals(m.fba().Z, flux)
        m = JsonModel.from_model(m).to_model()
        self.assertEquals(m.fba().Z, flux)
