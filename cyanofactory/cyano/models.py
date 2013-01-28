from django.db import models
import datetime
import math
from django.utils import timezone
from django.utils.safestring import mark_safe

from biosql.helpers import get_database_item
from db_xref.views import get_database_url_from_organism

class Gene:
    def __init__(self, organism, gene):
        """
        Creates a new Gene object.
        organism: Organism with the gene (Bioentry)
        gene: Name of the gene
        """
        self.organism = organism
        self.gene_name = gene
        self.record = get_database_item(organism)
        self.gene = None
        self.cds = None
        
        for item in self.record.features:
            if item.type == "gene":
                if "gene" in item.qualifiers:
                    if item.qualifiers["gene"][0] == self.gene_name:
                        self.gene = item
            elif item.type == "CDS":
                if "gene" in item.qualifiers:
                    if item.qualifiers["gene"][0] == self.gene_name:
                        self.cds = item
        
        
        
        #if self.gene == None:
            # error
        
        self.fieldsets = [#[
            #("Database", {'fields': [self.get_databases_html()]}),
            #("Moep", {"fields": ["aaa", "bbb"]})
            #]
            
            #('Type', ['model_type']),
            ('Name', [
                {'name': 'Name', 'data': gene},
                {'name': 'Synonyms', 'data': lambda : ", ".join(self.gene.qualifiers.get("gene_synonym", [""])[0].split("; "))},
                {'name': 'Cross References', 'data' : "", 'name': 'protein_complexes'},
                ]),
            ('Structure', [
                #{'name': "Structure", 'data': self.get_as_html_structure_local(100, 500)},
                {'name': 'Sequence', 'data': self.get_as_html_sequence()},
                {'name': 'Synonyms', 'data': lambda : ", ".join(self.gene.qualifiers.get("gene_synonym", [""])[0].split("; "))},
                {'name': 'Cross References', 'data' : "", 'name': 'protein_complexes'},
                ]),
            ('Comments', [
                {"name": "References",
                 "data": mark_safe(self.get_databases_html())}]),
            #('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]}),
            ]
    
    def get_databases_html(self):
        result = ""
        for x in self.gene.qualifiers["db_xref"]:
            result += '<a href="%s">%s</a>' % (get_database_url_from_organism(x), x)
            result += '<br />'
        return result

    def get_name(self):
        return self.get_name()
    
    def get_as_html_sequence(self):
        from cyano.helpers import format_sequence_as_html
        
        #direction = CHOICES_DIRECTION[[x[0] for x in CHOICES_DIRECTION].index(self.direction)][1]        
        
        #return 'Chromosome: <a href="%s">%s</a>, Coordinate: %s (nt), Length: %s (nt), Direction: %s, G/C content: %.1f%%, Sequence: %s' % (
        #    self.chromosome.get_absolute_url(), self.chromosome.wid, 
        #    self.coordinate, self.length, direction, 
        #    self.get_gc_content() * 100,
        #    format_sequence_as_html(self.get_sequence()))
        #return 'Sequence: %s' %        
        return 'Sequence: %s' % (format_sequence_as_html(self.gene.extract(self.record.seq)))

        
    class Meta:
        fieldsets = ["test", 
        #    ('Type', {'fields': ['model_type']}),
        #    ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), 
        #    ('Classification', {'fields': ['type']}), 
        #    ('Content', {'fields': [
        #        {'verbose_name': 'Metabolites (mM)', 'name': 'biomass_compositions'},
        #        {'verbose_name': 'Protein monomers', 'name': 'protein_monomers'},
        #        {'verbose_name': 'Protein complexes', 'name': 'protein_complexes'},
        #        ]}),
        #    ('Comments', {'fields': ['comments', 'references']}),
        #    ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]}),
        ]       
        field_list = []
        template = "cyano/detail.html"

class Genes():
    class Meta:
        field_list = [
            ["gene", "Gene"],
            ["locus_tag", "Locus Tag"],
            ["function", "Function"],
        ]
        field_url = "gene"
        template = "cyano/list.html"

# Create your models here.
#class Poll(models.Model):
#    question = models.CharField(max_length = 200)
#    pub_date = models.DateTimeField("data published")
#    
#    def was_published_recently(self):
#        return self.pub_date >= timezone.now() - datetime.timedelta(days = 1)
#    was_published_recently.admin_order_field = "pub_date"
#    was_published_recently.boolean = True
#    was_published_recently.short_description = "Published recently?"
#    
#    def __unicode__(self):
#        return self.question

#class Choice(models.Model):
#    poll = models.ForeignKey(Poll)
#    choice = models.CharField(max_length = 200)
#    votes = models.IntegerField()
#    
#    def __unicode__(self):
#        return self.choice
