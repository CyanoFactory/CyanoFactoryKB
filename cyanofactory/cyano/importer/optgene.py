from importer import Importer
from importer import ParserError
from django.core.exceptions import ObjectDoesNotExist
import cyano.models as models
from cyano.helpers import format_field_detail_view, render_queryset_to_response,\
    objectToQuerySet, format_field_detail_view_diff
from django.db.models.related import RelatedObject
from django.template.defaultfilters import capfirst
from bioparser.optgene import OptGeneParser
import re
import StringIO

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
