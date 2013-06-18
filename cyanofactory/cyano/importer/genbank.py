from importer import Importer
from importer import ParserError
from django.core.exceptions import ObjectDoesNotExist
import cyano.models as models
from cyano.helpers import format_field_detail_view, render_queryset_to_response,\
    objectToQuerySet, format_field_detail_view_diff
from django.db.models.related import RelatedObject
from django.template.defaultfilters import capfirst
import re
import StringIO

class GenbankImporter(Importer):
    def __init__(self):
        super(Importer, self).__init__()
        self.sequences = []

    def load(self, filename):
        self.filename = filename
        with open(filename) as f:
            self.data = list(for record in SeqIO.parse(f, "genbank"))

    def preview(self, request, species_wid):
        #super(Importer, self).preview(request, species_wid)
        
        try:
            species = models.Species.objects.get(wid = species_wid)
        except ObjectDoesNotExist:
            species = models.Species(wid = species_wid, name = request.POST["new_species"])
        
        fields = []
        
        for feature in self.data.features:
            qualifiers = feature.qualifiers
            wid = qualifiers["locus_tag"]
            name = qualifiers["gene"]
            comments = qualifiers["note"]
            
            try:
                old_gene = models.Gene.objects.get(wid = wid)
                new_gene = models.Gene.objects.get(wid = wid)
            except ObjectDoesNotExist:
                old_gene = None
                new_gene = models.Gene(wid = wid)     
            
            new_gene.name = name
            new_gene.comments = comments
            new_gene.wid = wid
            
            if old_gene == None:
                upd_message = "Create Gene"
            else:
                upd_message = "Update Gene"

            fieldsets = [
                (upd_message + " " + wid, {'fields': [
                    {'verbose_name': 'Name', 'name': 'name'},
                    {'verbose_name': 'Comments', 'name': 'comments'}]})
                ]
        
            #filter out type, metadata
            fieldset_names = [x[0] for x in fieldsets]
            new_gene.model_type = models.TableMeta.objects.get(name = "Gene")
            if 'Type' in fieldset_names:
                idx = fieldset_names.index('Type')
                del fieldsets[idx]
        
            model = models.Gene
            #filter out empty fields
            rmfieldsets = []
            for idx in range(len(fieldsets)):
                rmfields = []
                for idx2 in range(len(fieldsets[idx][1]['fields'])):
                    if isinstance(fieldsets[idx][1]['fields'][idx2], dict):
                        field_name = fieldsets[idx][1]['fields'][idx2]['name']
                        verbose_name = fieldsets[idx][1]['fields'][idx2]['verbose_name']
                    else:
                        field_name = fieldsets[idx][1]['fields'][idx2]
                        field = model._meta.get_field_by_name(field_name)[0]
                        if isinstance(field, RelatedObject):
                            verbose_name = capfirst(field.get_accessor_name())
                            pass
                        else:
                            verbose_name = field.verbose_name
                        
                    old_data, new_data = format_field_detail_view_diff(species, old_gene, new_gene, field_name, request.user.is_anonymous())
                    
                    if not old_data and not new_data:
                        rmfields = [idx2] + rmfields
                    
                    fieldsets[idx][1]['fields'][idx2] = {'verbose_name': verbose_name.replace(" ", '&nbsp;').replace("-", "&#8209;"), 'new_data': new_data, 'old_data': old_data}
                for idx2 in rmfields:
                    del fieldsets[idx][1]['fields'][idx2]
                if len(fieldsets[idx][1]['fields']) == 0:
                    rmfieldsets = [idx] + rmfieldsets
            for idx in rmfieldsets:
                del fieldsets[idx]
                
            fields.append(fieldsets[0])
        
        return render_queryset_to_response(
            species = species,
            request = request, 
            models = [model],
            queryset = objectToQuerySet(new_gene),
            template = 'cyano/preview.html', 
            data = {
                'model_type': "Chromosome",
                'model': model,
                'fieldsets': fields,
                'message': request.GET.get('message', ''),
                'filename': self.filename,
                'species_wid': species_wid
                })
    
    def submit(self, request, species_wid):
        revdetail = models.RevisionDetail()
        revdetail.user = request.user.profile
        revdetail.reason = request.POST["reason"]
        
        try:
            species = models.Species.objects.get(wid = species_wid)
        except ObjectDoesNotExist:
            species = models.Species(wid = species_wid, name = request.POST["new_species"])
            species.save(revision_detail = revdetail)
            
        try:
            chromosome = models.Chromosome.objects.get(wid = "CHROMOSOME-1")
        except ObjectDoesNotExist:
            chromosome = models.Chromosome(wid = "CHROMOSOME-1")
            chromosome.save(revision_detail = revdetail)
        
        chromosome.species.add(species)
        chromosome.save(revision_detail = revdetail)
        
        for feature in self.data.features:
            qualifiers = feature.qualifiers
            wid = qualifiers["locus_tag"]
            name = qualifiers["gene"]
            comments = qualifiers["note"]
            
            try:
                gene = models.Gene.objects.get(wid = wid)
            except ObjectDoesNotExist:
                gene = models.Gene(wid = wid)
                gene.chromosome = chromosome
            
            gene.name = name
            gene.direction = 'f' if feature.location.start < feature.location.end else 'r'
            gene.coordinate = feature.location.start if gene.direction == 'f' else feature.location.end
            gene.length = abs(feature.location.start - feature.location.end)
            gene.comments = comments
            gene.wid = wid
            
            gene.save(revision_detail = revdetail)
            gene.species.add(species)
            gene.save(revision_detail = revdetail)
        
        return render_queryset_to_response(
            species = species,
            request = request,
            template = 'cyano/importDataResult.html', 
            data = {
                    'success': 'success',
                    'message': "Data imported",
                   })