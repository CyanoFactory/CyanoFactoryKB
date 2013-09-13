from cyano.tests import CyanoBaseTest

class DBXrefTest(CyanoBaseTest):
    valid_test_url = "/dbxref/EC:1.2.3.4"
    invalid_test_url = "/dbxref/INVALID:ITEM"
    
    def test_main_page_works(self):                
        with self.assertTemplateUsed("db_xref/index.html"):
            self.assertOK("/dbxref/")
    
    def test_redirect(self):
        self.assertRedirect("/dbxref/EC:1.2.3.4")
        self.assertRedirect("/dbxref/eC:1.2.3.4")
        self.assertRedirect("/dbxref/Ec:1.2.3.4")
        self.assertRedirect("/dbxref/ec:1.2.3.4")
    
    def test_error_page_invalid_database(self):
        with self.assertTemplateUsed("db_xref/error.html"):
            self.assertBadRequest(self.invalid_test_url)
    
    def test_mimetype_json(self):
        with self.assertTemplateUsed("db_xref/output.json"):
            self.assertMimeTypeJson(self.valid_test_url + "?format=json")
    
    def test_mimetype_xml(self):
        with self.assertTemplateUsed("db_xref/output.xml"):
            self.assertMimeTypeXml(self.valid_test_url + "?format=xml")
    
    def test_mimetype_text(self):
        with self.assertTemplateUsed("db_xref/output.txt"):
            self.assertMimeTypeText(self.valid_test_url + "?format=txt")
    
    def test_mimetype_html(self):
        with self.assertTemplateUsed("db_xref/output.html"):
            self.GET(self.valid_test_url + "?format=html")
    
    def test_mimetype_json_invalid(self):
        with self.assertTemplateUsed("db_xref/output.json"):
            self.assertMimeTypeJson(self.invalid_test_url + "?format=json", status = 400)
    
    def test_mimetype_xml_fail_invalid(self):
        with self.assertTemplateUsed("db_xref/output.xml"):
            self.assertMimeTypeXml(self.invalid_test_url + "?format=xml", status = 400)
    
    def test_mimetype_text_fail_invalid(self):
        with self.assertTemplateUsed("db_xref/output.txt"):
            self.assertMimeTypeText(self.invalid_test_url + "?format=txt", status = 400)
    
    def test_mimetype_html_fail_invalid(self):
        with self.assertTemplateUsed("db_xref/output.html"):
            self.assertBadRequest(self.invalid_test_url + "?format=html")

    def test_error_invalid_format(self):
        with self.assertTemplateUsed("db_xref/error.html"):
            self.assertBadRequest(self.valid_test_url + "?format=invalid")
            self.assertBadRequest(self.invalid_test_url + "?format=invalid")
    
    def test_invalid_urls(self):
        self.assertNotFound("db_xref/INVALID:")
        self.assertNotFound("db_xref/EC:")
        self.assertNotFound("db_xref/:ABC")
