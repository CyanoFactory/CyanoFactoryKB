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

    def assertOK(self, url):
        """Fail if result code is not 200"""
        self.GET(url, status = 200)

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
        from management.commands.autocreateinitial import Command as autocreate
        autocreate().handle()
        
        # Create meta data tables
        from management.commands.create_meta import Command as create_meta
        create_meta().handle()
        
        # Create an entry
        self.setUpSpecies()
    
    def setUpSpecies(self):
        rev_detail = cmodels.RevisionDetail(
            user = cmodels.UserProfile.objects.get(user__username = "admin"),
            reason = "Unit Testing")
        rev_detail.save()

        species = cmodels.Species(wid = "UnitTest-Species")
        species.name = "Species for Unit Test"
        species.genetic_code = '11'
        species.save(rev_detail)
        
        chromosome = cmodels.Chromosome(wid = "UnitTest-Chromosome")
        chromosome.name = "Chromosome for Unit Test"
        chromosome.sequence = "AAAATTTTCCCCGGGGATCGATCGATCG"
        chromosome.length = len(chromosome.sequence)
        chromosome.save(rev_detail)
        chromosome.species.add(species)

        gene = cmodels.Gene(wid = "Test-Gene")
        gene.name = "Gene for Unit Test"
        gene.symbol = "xyzA"
        gene.direction = 'f'
        gene.length = 10
        gene.coordinate = 5
        gene.chromosome = chromosome
        gene.save(rev_detail)
        gene.species.add(species)
    
    def test_page_works(self):
        with self.assertTemplateUsed("cyano/index.html"): 
            self.assertOK("/")
    
    def test_permissions_species(self):
        self.assertForbidden("/UnitTest-Species")
        self.assertForbidden("/UnitTest-Species/")
        
        self.client.login(username = "gabriel", password = "aaa")
        self.assertForbidden("/UnitTest-Species")
        self.assertForbidden("/UnitTest-Species/")
        
        permission = cmodels.Permission.get_by_name("READ_NORMAL")
        
        user_perm = cmodels.UserPermission()
        user_perm.user = cmodels.UserProfile.objects.get(user__username = "gabriel")
        user_perm.entry = cmodels.Species.objects.get(wid = "UnitTest-Species")
        user_perm.save()
        user_perm.allow.add(permission)
        
        with self.assertTemplateUsed("cyano/species.html"): 
            self.assertOK("/UnitTest-Species")
            self.assertOK("/UnitTest-Species/")
        
        user_perm.deny.add(permission)
        
        self.assertForbidden("/UnitTest-Species/")
    
    def test_species_not_found(self):
        with self.assertTemplateUsed("cyano/error.html"):
            self.assertNotFound("/Wrong-Species")
    
    def test_logic_works(self):
        self.assertFalse(True == False, "Hardware error")
