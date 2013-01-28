# This is an auto-generated Django model module.
# You'll have to do the following manually to clean this up:
#     * Rearrange models' order
#     * Make sure each model has one field with primary_key=True
# Feel free to rename the models, but don't rename db_table values or field names.
#
# Also note: You'll have to insert the output of 'django-admin.py sqlcustom [appname]'
# into your database.

from django.db import models
from django.core.urlresolvers import reverse

class Ontology(models.Model):
    ontology_id = models.IntegerField(primary_key=True)
    name = models.CharField(max_length=32, unique=True)
    definition = models.TextField(blank=True)
    class Meta:
        db_table = u'ontology'

class Term(models.Model):
    term_id = models.IntegerField(primary_key=True)
    name = models.CharField(max_length=255)
    definition = models.TextField(blank=True)
    identifier = models.CharField(max_length=40, unique=True, blank=True)
    is_obsolete = models.CharField(max_length=1, blank=True)
    ontology = models.ForeignKey(Ontology)
    class Meta:
        db_table = u'term'

class Dbxref(models.Model):
    dbxref_id = models.IntegerField(primary_key=True)
    dbname = models.CharField(max_length=40)
    accession = models.CharField(max_length=128)
    version = models.IntegerField()
    class Meta:
        db_table = u'dbxref'

class Biodatabase(models.Model):
    biodatabase_id = models.IntegerField(primary_key=True)
    name = models.CharField(max_length=128, unique=True)
    authority = models.CharField(max_length=128, blank=True)
    description = models.TextField(blank=True)
    class Meta:
        db_table = u'biodatabase'

class Taxon(models.Model):
    taxon_id = models.IntegerField(primary_key=True)
    ncbi_taxon_id = models.IntegerField(unique=True, null=True, blank=True)
    parent_taxon_id = models.IntegerField(null=True, blank=True)
    node_rank = models.CharField(max_length=32, blank=True)
    genetic_code = models.SmallIntegerField(null=True, blank=True)
    mito_genetic_code = models.SmallIntegerField(null=True, blank=True)
    left_value = models.IntegerField(unique=True, null=True, blank=True)
    right_value = models.IntegerField(unique=True, null=True, blank=True)
    class Meta:
        db_table = u'taxon'

class Bioentry(models.Model):
    bioentry_id = models.IntegerField(primary_key=True)
    biodatabase = models.ForeignKey(Biodatabase)
    taxon = models.ForeignKey(Taxon, null=True, blank=True)
    name = models.CharField(max_length=40)
    accession = models.CharField(max_length=128)
    identifier = models.CharField(max_length=40, blank=True)
    division = models.CharField(max_length=6, blank=True)
    description = models.TextField(blank=True)
    version = models.IntegerField()
    
    def get_scientific_name(self):
        return TaxonName.objects.get(taxon_id = self.taxon_id, name_class = "scientific name").name
    
    def get_absolute_url(self):
        return reverse('cyano.views.organism', args=[self.name])
    
    def get_genes_url(self):
        return reverse('cyano.views.genes', args=[self.name])

    def get_gene_url(self, gene):
        return reverse('cyano.views.gene', kwargs={'organism': self.name,
                                                    'gene': gene})

    class Meta:
        db_table = u'bioentry'

class Seqfeature(models.Model):
    seqfeature_id = models.IntegerField(primary_key=True)
    bioentry = models.ForeignKey(Bioentry, related_name='seqfeature_bioentries')
    type_term = models.ForeignKey(Term, related_name='seqfeature_type_terms')
    source_term = models.ForeignKey(Term, related_name='seqfeature_source_terms')
    display_name = models.CharField(max_length=64, blank=True)
    rank = models.IntegerField()
    class Meta:
        db_table = u'seqfeature'
    
    @staticmethod
    def filter_by_organism(organism):
        return Seqfeature.objects.filter(bioentry__biodatabase__name = organism.biodatabase.name).filter(bioentry = organism)
        #             .filter(bioentry = organism)
        #Seqfeature.objects.filter(type_term__name = "gene").filter(bioentry__biodatabase__name = database_name)
        #return self.objects.filter(type_term__name = "gene")

class TermDbxref(models.Model):
    fake_id = models.IntegerField(primary_key=True)
    term = models.ForeignKey(Term)
    dbxref = models.ForeignKey(Dbxref)
    rank = models.IntegerField(null=True, blank=True)
    class Meta:
        db_table = u'term_dbxref'

class TermRelationship(models.Model):
    term_relationship_id = models.IntegerField(primary_key=True)
    subject_term = models.ForeignKey(Term, related_name='term_relationship_subject_terms')
    predicate_term = models.ForeignKey(Term, related_name='term_relationship_predicate_terms')
    object_term = models.ForeignKey(Term, related_name='term_relationship_object_terms')
    ontology = models.ForeignKey(Ontology)
    class Meta:
        db_table = u'term_relationship'

class BioentryDbxref(models.Model):
    fake_id = models.IntegerField(primary_key=True)
    bioentry = models.ForeignKey(Bioentry)
    dbxref = models.ForeignKey(Dbxref)
    rank = models.IntegerField(null=True, blank=True)
    class Meta:
        db_table = u'bioentry_dbxref'

class Biosequence(models.Model):
    bioentry = models.ForeignKey(Bioentry, primary_key=True)
    version = models.IntegerField(null=True, blank=True)
    length = models.IntegerField(null=True, blank=True)
    alphabet = models.CharField(max_length=10, blank=True)
    seq = models.TextField(blank=True)
    class Meta:
        db_table = u'biosequence'

class BioentryQualifierValue(models.Model):
    fake_id = models.IntegerField(primary_key=True)
    bioentry = models.ForeignKey(Bioentry)
    term = models.ForeignKey(Term)
    value = models.TextField(blank=True)
    rank = models.IntegerField()
    class Meta:
        db_table = u'bioentry_qualifier_value'

class Location(models.Model):
    location_id = models.IntegerField(primary_key=True)
    seqfeature = models.ForeignKey(Seqfeature)
    dbxref = models.ForeignKey(Dbxref, null=True, blank=True)
    term = models.ForeignKey(Term, null=True, blank=True)
    start_pos = models.IntegerField(null=True, blank=True)
    end_pos = models.IntegerField(null=True, blank=True)
    strand = models.IntegerField()
    rank = models.IntegerField()
    class Meta:
        db_table = u'location'

class Reference(models.Model):
    reference_id = models.IntegerField(primary_key=True)
    dbxref = models.ForeignKey(Dbxref, unique=True, null=True, blank=True)
    location = models.TextField()
    title = models.TextField(blank=True)
    authors = models.TextField(blank=True)
    crc = models.CharField(max_length=32, unique=True, blank=True)
    class Meta:
        db_table = u'reference'

class DbxrefQualifierValue(models.Model):
    fake_id = models.IntegerField(primary_key=True)
    dbxref = models.ForeignKey(Dbxref)
    term = models.ForeignKey(Term)
    rank = models.IntegerField()
    value = models.TextField(blank=True)
    class Meta:
        db_table = u'dbxref_qualifier_value'

class LocationQualifierValue(models.Model):
    fake_id = models.IntegerField(primary_key=True)
    location = models.ForeignKey(Location)
    term = models.ForeignKey(Term)
    value = models.CharField(max_length=255)
    int_value = models.IntegerField(null=True, blank=True)
    class Meta:
        db_table = u'location_qualifier_value'

class SeqfeatureDbxref(models.Model):
    fake_id = models.IntegerField(primary_key=True)
    seqfeature = models.ForeignKey(Seqfeature)
    dbxref = models.ForeignKey(Dbxref)
    rank = models.IntegerField(null=True, blank=True)
    class Meta:
        db_table = u'seqfeature_dbxref'

class SeqfeaturePath(models.Model):
    fake_id = models.IntegerField(primary_key=True)
    object_seqfeature = models.ForeignKey(Seqfeature, related_name='seqfeature_object_seqfeatures')
    subject_seqfeature = models.ForeignKey(Seqfeature, related_name='seqfeature_subject_seqfeatures')
    term = models.ForeignKey(Term, related_name='seqfeature_path_terms')
    distance = models.IntegerField(null=True, blank=True)
    class Meta:
        db_table = u'seqfeature_path'

class SeqfeatureQualifierValue(models.Model):
    fake_id = models.IntegerField(primary_key=True)
    seqfeature = models.ForeignKey(Seqfeature)
    term = models.ForeignKey(Term)
    rank = models.IntegerField()
    value = models.TextField()
    class Meta:
        db_table = u'seqfeature_qualifier_value'

class SeqfeatureRelationship(models.Model):
    seqfeature_relationship_id = models.IntegerField(primary_key=True)
    object_seqfeature = models.ForeignKey(Seqfeature, related_name='seqfeature_relationship_object_seqfeatures')
    subject_seqfeature = models.ForeignKey(Seqfeature, related_name='seqfeature_relationship_subject_subfeatures')
    term = models.ForeignKey(Term, related_name='seqfeature_relationship_terms')
    rank = models.IntegerField(null=True, blank=True)
    class Meta:
        db_table = u'seqfeature_relationship'

class TaxonName(models.Model):
    fake_id = models.IntegerField(primary_key=True)
    taxon = models.ForeignKey(Taxon)
    name = models.CharField(max_length=255)
    name_class = models.CharField(max_length=32)
    class Meta:
        db_table = u'taxon_name'

class TermPath(models.Model):
    term_path_id = models.IntegerField(primary_key=True)
    subject_term = models.ForeignKey(Term, related_name='term_path_subject_terms')
    predicate_term = models.ForeignKey(Term, related_name='term_path_predicate_terms')
    object_term = models.ForeignKey(Term, related_name='term_path_object_terms')
    ontology = models.ForeignKey(Ontology, related_name='term_path_ontologies')
    distance = models.IntegerField(null=True, blank=True)
    class Meta:
        db_table = u'term_path'

class TermRelationshipTerm(models.Model):
    term_relationship = models.ForeignKey(TermRelationship, primary_key=True)
    term = models.ForeignKey(Term, unique=True)
    class Meta:
        db_table = u'term_relationship_term'

class BioentryPath(models.Model):
    fake_id = models.IntegerField(primary_key=True)
    object_bioentry = models.ForeignKey(Bioentry, related_name='bioentry_path_object_bioentries')
    subject_bioentry = models.ForeignKey(Bioentry, related_name='bioentry_path_subject_bioentries')
    term = models.ForeignKey(Term, related_name='bioentry_path_terms')
    distance = models.IntegerField(null=True, blank=True)
    class Meta:
        db_table = u'bioentry_path'

class TermSynonym(models.Model):
    fake_id = models.IntegerField(primary_key=True)
    synonym = models.CharField(max_length=255)
    term = models.ForeignKey(Term)
    class Meta:
        db_table = u'term_synonym'

class BioentryReference(models.Model):
    fake_id = models.IntegerField(primary_key=True)
    bioentry = models.ForeignKey(Bioentry)
    reference = models.ForeignKey(Reference)
    start_pos = models.IntegerField(null=True, blank=True)
    end_pos = models.IntegerField(null=True, blank=True)
    rank = models.IntegerField()
    class Meta:
        db_table = u'bioentry_reference'

class BioentryRelationship(models.Model):
    bioentry_relationship_id = models.IntegerField(primary_key=True)
    object_bioentry = models.ForeignKey(Bioentry, related_name='bioentry_relationship_object_bioentries')
    subject_bioentry = models.ForeignKey(Bioentry, related_name='bioentry_relationship_subject_bioentries')
    term = models.ForeignKey(Term, related_name='bioentry_relationship_terms')
    rank = models.IntegerField(null=True, blank=True)
    class Meta:
        db_table = u'bioentry_relationship'

class Comment(models.Model):
    comment_id = models.IntegerField(primary_key=True)
    bioentry = models.ForeignKey(Bioentry)
    comment_text = models.TextField()
    rank = models.IntegerField()
    class Meta:
        db_table = u'comment'

