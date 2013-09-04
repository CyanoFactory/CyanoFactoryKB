from importer import Importer
from bioparser.optgene import OptGeneParser

class OptGeneImporter(Importer):
    def __init__(self):
        super(Importer, self).__init__()

    def load(self, filename):
        self.filename = filename
        self.data = OptGeneParser(filename)

    def preview(self, request, species_wid):
        pass
        #super(Importer, self).preview(request, species_wid)
    
    def submit(self, request, species_wid):
        pass
