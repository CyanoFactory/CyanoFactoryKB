from cyano.tests import CyanoBaseTest
from django.core import management

def setup():
    import sys

    # Create users and groups
    management.call_command('autocreateinitial')
    # Import boehringer (reads from stdin)
    sys.stdin = open("boehringer/fixtures/data.dump")
    management.call_command('restoredb')
    sys.stdin.close()
    sys.stdin = sys.__stdin__    

class BoehringerStandaloneGuestTest(CyanoBaseTest):
    """Tests standalone feature under /boehringer/ as guest"""
    
    def test_main_page_guest(self):
        """Main page as guest"""
        self.assertLoginRequired("boehringer/")

    def test_legacy_page_guest(self):
        """Main page as guest"""
        self.assertLoginRequired("boehringer/legacy/")

class BoehringerStandaloneTest(CyanoBaseTest):
    """Tests standalone feature under /boehringer/"""
    def setUp(self):
        self.doLogin()
    
    def test_main_page_works(self):
        """Main page loads"""        
        with self.assertTemplateUsed("boehringer/index.html"):
            self.assertOK("boehringer/")
    
    def test_legacy_page_works(self):
        """Legacy page loads"""
        with self.assertTemplateUsed("boehringer/legacy.html"):
            self.assertOK("boehringer/legacy/")
    
    def test_default_search(self):
        """default search results are correct"""
        with self.assertTemplateUsed("boehringer/index.html"):
            resp = self.assertOK("boehringer/")
            self.assertContains(resp, "4 hits - 1 miss")
            self.assertContains(resp, "1 hit - 0 misses")
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
