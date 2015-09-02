# This is an auto-generated Django model module.
# You'll have to do the following manually to clean this up:
#     * Rearrange models' order
#     * Make sure each model has one field with primary_key=True
# Feel free to rename the models, but don't rename db_table values or field names.
#
# Also note: You'll have to insert the output of 'django-admin.py sqlcustom [appname]'
# into your database.
from __future__ import unicode_literals

from django.db import models
from django.db.models import options
import dbarray

options.DEFAULT_NAMES = options.DEFAULT_NAMES + ('db_interaction_div',)


class Abstracts(models.Model):
    abstract_id = models.CharField(max_length=250)
    publication_date = models.SmallIntegerField(null=True, blank=True)
    publication_source = models.TextField(blank=True)
    linkout_url = models.TextField(blank=True)
    title = models.TextField()
    body = models.TextField()

    class Meta:
        db_interaction_div = "string"
        db_table = 'abstracts'


class Actions(models.Model):
    item_id_a = models.IntegerField()
    item_id_b = models.IntegerField()
    mode = models.CharField(max_length=255)
    action = models.CharField(max_length=255)
    a_is_acting = models.BooleanField()
    score = models.SmallIntegerField()

    class Meta:
        db_interaction_div = "string"
        db_table = 'actions'


class ActionsSets(models.Model):
    item_id_a = models.IntegerField()
    item_id_b = models.IntegerField()
    mode = models.CharField(max_length=255)
    sources = models.TextField()

    class Meta:
        db_interaction_div = "string"
        db_table = 'actions_sets'


class BestCombinedScoresOrthgroups(models.Model):
    orthgroup_id = models.IntegerField(primary_key=True)
    best_score = models.IntegerField()

    class Meta:
        db_interaction_div = "string"
        db_table = 'best_combined_scores_orthgroups'


class BestCombinedScoresProteins(models.Model):
    protein_id = models.IntegerField(primary_key=True)
    best_score = models.IntegerField()

    class Meta:
        db_interaction_div = "string"
        db_table = 'best_combined_scores_proteins'


class Collections(models.Model):
    collection_id = models.IntegerField(primary_key=True)
    pubmed_id = models.IntegerField()
    comment = models.CharField(max_length=1000, blank=True)

    class Meta:
        db_interaction_div = "string"
        db_table = 'collections'


class EvidenceTransfers(models.Model):
    target_protein_id_a = models.IntegerField()
    target_protein_id_b = models.IntegerField()
    source_protein_id_a = models.IntegerField()
    source_protein_id_b = models.IntegerField()
    transfer_score_c1 = models.SmallIntegerField()
    transfer_score_c2 = models.SmallIntegerField()
    source_ascore = models.SmallIntegerField()
    source_escore = models.SmallIntegerField()
    source_dscore = models.SmallIntegerField()
    source_tscore = models.SmallIntegerField()

    class Meta:
        db_interaction_div = "string"
        db_table = 'evidence_transfers'


class Funccats(models.Model):
    funccat_id = models.CharField(max_length=5)
    funccat_description = models.CharField(max_length=100)

    class Meta:
        db_interaction_div = "string"
        db_table = 'funccats'


class FusionEvidence(models.Model):
    target_protein_id_a = models.IntegerField()
    target_protein_id_b = models.IntegerField()
    source_protein = models.IntegerField()
    source_species = models.IntegerField()
    transfer_score_c1 = models.SmallIntegerField()
    transfer_score_c2 = models.SmallIntegerField()
    fusion_score = models.SmallIntegerField()

    class Meta:
        db_interaction_div = "string"
        db_table = 'fusion_evidence'


class Genes(models.Model):
    gene_id = models.IntegerField(primary_key=True)
    gene_external_id = models.CharField(max_length=100)
    start_position_on_contig = models.IntegerField()
    end_position_on_contig = models.IntegerField()
    protein_size = models.IntegerField()

    class Meta:
        db_interaction_div = "string"
        db_table = 'genes'


class GenesProteins(models.Model):
    protein_id = models.IntegerField()
    gene_id = models.IntegerField()

    class Meta:
        db_interaction_div = "string"
        db_table = 'genes_proteins'


class ItemsAbstracts(models.Model):
    protein_id = models.IntegerField()
    abstract_id = models.CharField(max_length=20)
    name = models.CharField(max_length=100)
    abstract_length = models.IntegerField()
    mesh_id = models.IntegerField()

    class Meta:
        db_interaction_div = "string"
        db_table = 'items_abstracts'


class Meshterms(models.Model):
    mesh_id = models.IntegerField(primary_key=True)
    description = models.CharField(max_length=255)

    class Meta:
        db_interaction_div = "string"
        db_table = 'meshterms'


class Orthgroups(models.Model):
    orthgroup_id = models.IntegerField(primary_key=True)
    orthgroup_external_id = models.CharField(max_length=20)
    description = models.CharField(max_length=1000)
    protein_count = models.IntegerField()
    species_count = models.SmallIntegerField()

    class Meta:
        db_interaction_div = "string"
        db_table = 'orthgroups'


class OrthgroupsAbstracts(models.Model):
    orthgroup_id = models.IntegerField()
    abstract_id = models.CharField(max_length=250)

    class Meta:
        db_interaction_div = "string"
        db_table = 'orthgroups_abstracts'


class OrthgroupsFunccats(models.Model):
    orthgroup_id = models.IntegerField()
    funccat_id = models.IntegerField()

    class Meta:
        db_interaction_div = "string"
        db_table = 'orthgroups_funccats'


class OrthgroupsSets(models.Model):
    orthgroup_id = models.IntegerField()
    set_id = models.CharField(max_length=60)
    is_database_set = models.BooleanField()

    class Meta:
        db_interaction_div = "string"
        db_table = 'orthgroups_sets'


class OrthgroupsSpecies(models.Model):
    orthgroup_id = models.IntegerField()
    species_id = models.IntegerField()
    count = models.IntegerField()

    class Meta:
        db_interaction_div = "string"
        db_table = 'orthgroups_species'


class ProteinImageMatch(models.Model):
    protein_id = models.IntegerField()
    image_id = models.CharField(max_length=250)
    identity = models.DecimalField(null=True, max_digits=65535, decimal_places=65535, blank=True)
    source = models.CharField(max_length=10)
    start_position_on_protein = models.IntegerField(null=True, blank=True)
    end_position_on_protein = models.IntegerField(null=True, blank=True)
    annotation = models.CharField(max_length=50, blank=True)

    class Meta:
        db_interaction_div = "string"
        db_table = 'protein_image_match'


class NodeNodeLinks(models.Model):
    node_id_a = models.IntegerField()
    node_type_b = models.IntegerField(null=True, blank=True)
    node_id_b = models.IntegerField()
    combined_score = models.SmallIntegerField(null=True, blank=True)
    evidence_scores = dbarray.TextArrayField()
    evidence = {}

    class Meta:
        db_interaction_div = "string"
        db_table = 'node_node_links'


class Proteins(models.Model):
    protein_id = models.IntegerField(primary_key=True)
    protein_external_id = models.CharField(max_length=100)
    species_id = models.IntegerField()
    protein_checksum = models.CharField(max_length=16)
    protein_size = models.IntegerField()
    annotation = models.CharField(max_length=600)
    preferred_name = models.CharField(max_length=50)
    annotation_word_vectors = models.TextField(blank=True)  # This field type is a guess.

    class Meta:
        db_interaction_div = "string"
        db_table = 'proteins'


class ProteinsMeshterms(models.Model):
    mesh_id = models.IntegerField()
    protein_id = models.IntegerField()

    class Meta:
        db_interaction_div = "string"
        db_table = 'proteins_meshterms'


class ProteinsNames(models.Model):
    protein_name = models.CharField(max_length=60)
    protein_id = models.IntegerField(primary_key=True)
    species_id = models.IntegerField()
    source = models.CharField(max_length=100)
    is_preferred_name = models.NullBooleanField(null=True, blank=True)
    linkout = models.CharField(max_length=15, blank=True)

    class Meta:
        db_interaction_div = "string"
        db_table = 'proteins_names'


class ProteinsOrthgroups(models.Model):
    orthgroup_id = models.IntegerField()
    protein_id = models.IntegerField()
    protein_external_id = models.CharField(max_length=50)
    species_id = models.IntegerField()
    start_position = models.IntegerField()
    end_position = models.IntegerField()
    preferred_name = models.CharField(max_length=50)
    protein_annotation = models.CharField(max_length=100)
    preferred_linkout_url = models.CharField(max_length=150, blank=True)

    class Meta:
        db_interaction_div = "string"
        db_table = 'proteins_orthgroups'


class ProteinsSequences(models.Model):
    # protein_id = models.IntegerField()
    protein_id = models.IntegerField(primary_key=True)
    sequence = models.TextField()

    class Meta:
        db_interaction_div = "string"
        db_table = 'proteins_sequences'


class ProteinsSmartlinkouts(models.Model):
    protein_id = models.IntegerField()
    protein_size = models.IntegerField()
    smart_url = models.CharField(max_length=2000)

    class Meta:
        db_interaction_div = "string"
        db_table = 'proteins_smartlinkouts'


class Runs(models.Model):
    run_id = models.IntegerField()
    species_id = models.IntegerField()
    contig_id = models.CharField(max_length=50)

    class Meta:
        db_interaction_div = "string"
        db_table = 'runs'


class RunsGenesProteins(models.Model):
    run_id = models.IntegerField()
    gene_id = models.IntegerField()
    protein_id = models.IntegerField()
    start_position_on_contig = models.IntegerField()
    end_position_on_contig = models.IntegerField()
    preferred_name = models.CharField(max_length=50)
    annotation = models.CharField(max_length=100)

    class Meta:
        db_interaction_div = "string"
        db_table = 'runs_genes_proteins'


class RunsOrthgroups(models.Model):
    run_id = models.IntegerField()
    orthgroup_id = models.IntegerField()

    class Meta:
        db_interaction_div = "string"
        db_table = 'runs_orthgroups'


class ScoreTypes(models.Model):
    score_id = models.SmallIntegerField()
    score_type = models.CharField(max_length=35)

    class Meta:
        db_interaction_div = "string"
        db_table = 'score_types'


class Sets(models.Model):
    set_id = models.CharField(max_length=60)
    collection_id = models.IntegerField()
    title = models.CharField(max_length=100, blank=True)
    comment = models.CharField(max_length=255, blank=True)
    url = models.CharField(max_length=255, blank=True)

    class Meta:
        db_interaction_div = "string"
        db_table = 'sets'


class SetsItems(models.Model):
    set_id = models.IntegerField()
    protein_id = models.IntegerField()
    species_id = models.IntegerField()
    is_database_set = models.BooleanField()

    class Meta:
        db_interaction_div = "string"
        db_table = 'sets_items'


class SetsPubmedrefs(models.Model):
    set_id = models.CharField(max_length=60)
    pubmed_id = models.CharField(max_length=20, blank=True)

    class Meta:
        db_interaction_div = "string"
        db_table = 'sets_pubmedrefs'


class Species(models.Model):
    species_id = models.IntegerField(primary_key=True)
    official_name = models.CharField(max_length=100)
    compact_name = models.CharField(max_length=100)
    kingdom = models.CharField(max_length=15)
    type = models.CharField(max_length=10)

    class Meta:
        db_interaction_div = "string"
        db_table = 'species'


class SpeciesNames(models.Model):
    species_id = models.IntegerField()
    species_name = models.IntegerField()
    official_name = models.CharField(max_length=255, blank=True)
    is_string_species = models.NullBooleanField(null=True, blank=True)

    class Meta:
        db_interaction_div = "string"
        db_table = 'species_names'


class SpeciesNodes(models.Model):
    species_id = models.IntegerField()
    species_name = models.CharField(max_length=255)
    position = models.IntegerField()
    size = models.IntegerField()

    class Meta:
        db_interaction_div = "string"
        db_table = 'species_nodes'


class ChemicalAliases(models.Model):
    chemical = models.CharField(max_length=50)
    alias = models.TextField()
    source = models.TextField()

    class Meta:
        db_interaction_div = "stitch"
        db_table = 'chemical_aliases'


class ChemicalChemicalLinks(models.Model):
    chemical1 = models.CharField(max_length=50)
    chemical2 = models.CharField(max_length=50)
    textmining = models.IntegerField()

    class Meta:
        db_interaction_div = "stitch"
        db_table = 'chemical_chemical_links'


class ChemicalChemicalLinksDetailed(models.Model):
    chemical1 = models.CharField(max_length=50)
    chemical2 = models.CharField(max_length=50)
    similarity = models.IntegerField()
    experimental = models.IntegerField()
    database = models.IntegerField()
    textmining = models.IntegerField()
    combined_score = models.IntegerField()

    class Meta:
        db_interaction_div = "stitch"
        db_table = 'chemical_chemical_links_detailed'


class ChemicalSources(models.Model):
    flat_chemical = models.CharField(max_length=50)
    stereo_chemical = models.CharField(max_length=50)
    source_name = models.CharField(max_length=50)
    source_id = models.CharField(max_length=50)

    class Meta:
        db_interaction_div = "stitch"
        db_table = 'chemical_sources'


class Chemicals(models.Model):
    chemical_id = models.IntegerField()
    preferred_name = models.TextField()
    short_name = models.TextField()
    mw = models.FloatField()
    drug = models.BooleanField()

    class Meta:
        db_interaction_div = "stitch"
        db_table = 'chemicals'


class ChemicalsInchikeys(models.Model):
    flat_chemical_id = models.CharField(max_length=50)
    stereo_chemical_id = models.CharField(max_length=50)
    source_cid = models.IntegerField()
    inchikey = models.TextField()

    class Meta:
        db_interaction_div = "stitch"
        db_table = 'chemicals_inchikeys'


class ProteinChemicalLinks(models.Model):
    chemical = models.CharField(max_length=50)
    protein = models.CharField(max_length=50)
    combined_score = models.IntegerField()

    class Meta:
        db_interaction_div = "stitch"
        # db_table = 'protein_chemical_links'
        db_table = 'node_node_links'


class ProteinChemicalLinksDetailed(models.Model):
    chemical = models.CharField(max_length=50)
    protein = models.CharField(max_length=50)
    experimental = models.IntegerField()
    prediction = models.IntegerField()
    database = models.IntegerField()
    textmining = models.IntegerField()
    combined_score = models.IntegerField()

    class Meta:
        db_interaction_div = "stitch"
        db_table = 'protein_chemical_links_detailed'


class ProteinChemicalLinksTransfer(models.Model):
    chemical = models.CharField(max_length=50)
    protein = models.CharField(max_length=50)
    experimental_direct = models.IntegerField()
    experimental_transferred = models.IntegerField()
    prediction_direct = models.IntegerField()
    prediction_transferred = models.IntegerField()
    database_direct = models.IntegerField()
    database_transferred = models.IntegerField()
    textmining_direct = models.IntegerField()
    textmining_transferred = models.IntegerField()
    combined_score = models.IntegerField()

    class Meta:
        db_interaction_div = "stitch"
        db_table = 'protein_chemical_links_transfer'


class Stitchactions(models.Model):
    item_id_a = models.CharField(max_length=50)
    item_id_b = models.CharField(max_length=50)
    mode = models.CharField(max_length=50)
    action = models.CharField(max_length=50, blank=True)
    a_is_acting = models.IntegerField()
    score = models.IntegerField()

    class Meta:
        db_interaction_div = "stitch"
        db_table = 'stitchactions'


class StitchNodeNode(models.Model):
    node_id_a = models.IntegerField()
    node_type_b = models.IntegerField()
    node_id_b = models.IntegerField()
    combined_score = models.SmallIntegerField()
    evidence_scores = dbarray.TextArrayField()

    class Meta:
        db_interaction_div = "stitch"
        db_table = "node_node_links"
