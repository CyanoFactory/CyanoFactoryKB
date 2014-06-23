# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'BasketComponent'
        db.create_table(u'cyano_basketcomponent', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('basket', self.gf('cyano.history.HistoryForeignKey')(related_name=u'components', to=orm['cyano.Basket'])),
            ('component', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', to=orm['cyano.SpeciesComponent'])),
            ('species', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', to=orm['cyano.Species'])),
        ))
        db.send_create_signal(u'cyano', ['BasketComponent'])

        # Adding model 'Basket'
        db.create_table(u'cyano_basket', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('user', self.gf('cyano.history.HistoryForeignKey')(related_name=u'baskets', to=orm['cyano.UserProfile'])),
            ('name', self.gf('django.db.models.fields.CharField')(default=u'', max_length=255)),
        ))
        db.send_create_signal(u'cyano', ['Basket'])


    def backwards(self, orm):
        # Deleting model 'BasketComponent'
        db.delete_table(u'cyano_basketcomponent')

        # Deleting model 'Basket'
        db.delete_table(u'cyano_basket')


    models = {
        u'auth.group': {
            'Meta': {'object_name': 'Group'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '80'}),
            'permissions': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['auth.Permission']", 'symmetrical': 'False', 'blank': 'True'})
        },
        u'auth.permission': {
            'Meta': {'ordering': "(u'content_type__app_label', u'content_type__model', u'codename')", 'unique_together': "((u'content_type', u'codename'),)", 'object_name': 'Permission'},
            'codename': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'content_type': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['contenttypes.ContentType']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'auth.user': {
            'Meta': {'object_name': 'User'},
            'date_joined': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'email': ('django.db.models.fields.EmailField', [], {'max_length': '75', 'blank': 'True'}),
            'first_name': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'groups': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'related_name': "u'user_set'", 'blank': 'True', 'to': u"orm['auth.Group']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'is_active': ('django.db.models.fields.BooleanField', [], {'default': 'True'}),
            'is_staff': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'is_superuser': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'last_login': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'last_name': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'password': ('django.db.models.fields.CharField', [], {'max_length': '128'}),
            'user_permissions': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'related_name': "u'user_set'", 'blank': 'True', 'to': u"orm['auth.Permission']"}),
            'username': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '30'})
        },
        u'contenttypes.contenttype': {
            'Meta': {'ordering': "('name',)", 'unique_together': "(('app_label', 'model'),)", 'object_name': 'ContentType', 'db_table': "'django_content_type'"},
            'app_label': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'model': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100'})
        },
        u'cyano.basket': {
            'Meta': {'object_name': 'Basket'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255'}),
            'user': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'baskets'", 'to': u"orm['cyano.UserProfile']"})
        },
        u'cyano.basketcomponent': {
            'Meta': {'object_name': 'BasketComponent'},
            'basket': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'components'", 'to': u"orm['cyano.Basket']"}),
            'component': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.SpeciesComponent']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'species': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.Species']"})
        },
        u'cyano.bindingsite': {
            'Meta': {'ordering': "[u'coordinate', u'length']", 'object_name': 'BindingSite'},
            'coordinate': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'direction': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'length': ('django.db.models.fields.PositiveIntegerField', [], {})
        },
        u'cyano.biomasscomposition': {
            'Meta': {'ordering': "[u'-concentration']", 'object_name': 'BiomassComposition'},
            'compartment': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'biomass_compositions'", 'to': u"orm['cyano.Compartment']"}),
            'concentration': ('django.db.models.fields.FloatField', [], {}),
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'})
        },
        u'cyano.chromosome': {
            'Meta': {'object_name': 'Chromosome', '_ormbases': [u'cyano.Genome']},
            'parent_ptr_genome': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "u'child_ptr_chromosome'", 'unique': 'True', 'primary_key': 'True', 'to': u"orm['cyano.Genome']"})
        },
        u'cyano.chromosomefeature': {
            'Meta': {'object_name': 'ChromosomeFeature', '_ormbases': [u'cyano.SpeciesComponent']},
            'parent_ptr_species_component': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "u'child_ptr_chromosome_feature'", 'unique': 'True', 'primary_key': 'True', 'to': u"orm['cyano.SpeciesComponent']"})
        },
        u'cyano.codon': {
            'Meta': {'ordering': "[u'sequence']", 'object_name': 'Codon'},
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'sequence': ('django.db.models.fields.CharField', [], {'max_length': '3'})
        },
        u'cyano.coenzymeparticipant': {
            'Meta': {'object_name': 'CoenzymeParticipant'},
            'coefficient': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'compartment': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.Compartment']"}),
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'metabolite': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'coenzyme_participants'", 'to': u"orm['cyano.Metabolite']"})
        },
        u'cyano.compartment': {
            'Meta': {'object_name': 'Compartment', '_ormbases': [u'cyano.SpeciesComponent']},
            'parent_ptr_species_component': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "u'child_ptr_compartment'", 'unique': 'True', 'primary_key': 'True', 'to': u"orm['cyano.SpeciesComponent']"})
        },
        u'cyano.crossreference': {
            'Meta': {'ordering': "[u'xid']", 'unique_together': "((u'xid', u'source'),)", 'object_name': 'CrossReference'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'source': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'xid': ('django.db.models.fields.CharField', [], {'max_length': '255'})
        },
        u'cyano.disulfidebond': {
            'Meta': {'object_name': 'DisulfideBond'},
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'protein_monomer': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'disulfide_bonds'", 'to': u"orm['cyano.ProteinMonomer']"}),
            'residue_1': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'residue_2': ('django.db.models.fields.PositiveIntegerField', [], {})
        },
        u'cyano.dnafootprint': {
            'Meta': {'object_name': 'DNAFootprint'},
            'binding': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '10', 'blank': 'True'}),
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'length': ('django.db.models.fields.PositiveIntegerField', [], {'null': 'True', 'blank': 'True'}),
            'region': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '10', 'blank': 'True'})
        },
        u'cyano.entry': {
            'Meta': {'ordering': "[u'wid']", 'object_name': 'Entry'},
            'comments': ('django.db.models.fields.TextField', [], {'default': "u''", 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'model_type': ('cyano.history.HistoryForeignKey', [], {'to': u"orm['cyano.TableMeta']"}),
            'name': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255', 'blank': 'True'}),
            'synonyms': ('cyano.history.HistoryManyToManyField', [], {'blank': 'True', 'related_name': "u'entry'", 'null': 'True', 'symmetrical': 'False', 'to': u"orm['cyano.Synonym']"}),
            'wid': ('django.db.models.fields.SlugField', [], {'max_length': '150'})
        },
        u'cyano.entrybasictextdata': {
            'Meta': {'object_name': 'EntryBasicTextData'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'value': ('django.db.models.fields.TextField', [], {'default': "u''"})
        },
        u'cyano.entrybooleandata': {
            'Meta': {'ordering': "[u'value']", 'object_name': 'EntryBooleanData'},
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'value': ('django.db.models.fields.BooleanField', [], {})
        },
        u'cyano.entrychardata': {
            'Meta': {'ordering': "[u'value', u'units']", 'object_name': 'EntryCharData'},
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'units': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255', 'blank': 'True'}),
            'value': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255', 'blank': 'True'})
        },
        u'cyano.entryfloatdata': {
            'Meta': {'ordering': "[u'value', u'units']", 'object_name': 'EntryFloatData'},
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'units': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255', 'blank': 'True'}),
            'value': ('django.db.models.fields.FloatField', [], {})
        },
        u'cyano.entrypositivefloatdata': {
            'Meta': {'ordering': "[u'value', u'units']", 'object_name': 'EntryPositiveFloatData'},
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'units': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255', 'blank': 'True'}),
            'value': ('django.db.models.fields.FloatField', [], {})
        },
        u'cyano.entrytextdata': {
            'Meta': {'ordering': "[u'value', u'units']", 'object_name': 'EntryTextData'},
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'units': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255', 'blank': 'True'}),
            'value': ('django.db.models.fields.TextField', [], {'default': "u''", 'blank': 'True'})
        },
        u'cyano.enzymeparticipant': {
            'Meta': {'object_name': 'EnzymeParticipant'},
            'compartment': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.Compartment']"}),
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'protein': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'enzyme_participants'", 'to': u"orm['cyano.Protein']"})
        },
        u'cyano.evidence': {
            'Meta': {'ordering': "[u'value', u'units']", 'object_name': 'Evidence'},
            'comments': ('django.db.models.fields.TextField', [], {'default': "u''", 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'is_experimentally_constrained': ('django.db.models.fields.BooleanField', [], {}),
            'media': ('django.db.models.fields.CharField', [], {'max_length': '255', 'null': 'True', 'blank': 'True'}),
            'pH': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'references': ('cyano.history.HistoryForeignKey', [], {'blank': 'True', 'related_name': "u'evidence'", 'null': 'True', 'to': u"orm['cyano.PublicationReference']"}),
            'species': ('django.db.models.fields.CharField', [], {'max_length': '255', 'null': 'True', 'blank': 'True'}),
            'species_component': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'related_name': "u'+'", 'blank': 'True', 'to': u"orm['cyano.SpeciesComponent']"}),
            'temperature': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'units': ('django.db.models.fields.CharField', [], {'max_length': '255', 'null': 'True', 'blank': 'True'}),
            'value': ('django.db.models.fields.TextField', [], {'default': "u''", 'blank': 'True'})
        },
        u'cyano.featureposition': {
            'Meta': {'object_name': 'FeaturePosition'},
            'chromosome': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'features'", 'to': u"orm['cyano.Genome']"}),
            'chromosome_feature': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'positions'", 'to': u"orm['cyano.ChromosomeFeature']"}),
            'coordinate': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'direction': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'length': ('django.db.models.fields.PositiveIntegerField', [], {})
        },
        u'cyano.gene': {
            'Meta': {'object_name': 'Gene', '_ormbases': [u'cyano.Molecule']},
            'amino_acid': ('cyano.history.HistoryForeignKey', [], {'blank': 'True', 'related_name': "u'genes'", 'null': 'True', 'on_delete': 'models.SET_NULL', 'to': u"orm['cyano.Metabolite']"}),
            'chromosome': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'genes'", 'to': u"orm['cyano.Genome']"}),
            'codons': ('cyano.history.HistoryManyToManyField', [], {'blank': 'True', 'related_name': "u'genes'", 'null': 'True', 'symmetrical': 'False', 'to': u"orm['cyano.Codon']"}),
            'coordinate': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'direction': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'expression': ('cyano.history.HistoryForeignKey', [], {'blank': 'True', 'related_name': "u'+'", 'null': 'True', 'to': u"orm['cyano.EntryPositiveFloatData']"}),
            'half_life': ('cyano.history.HistoryForeignKey', [], {'blank': 'True', 'related_name': "u'+'", 'null': 'True', 'to': u"orm['cyano.EntryPositiveFloatData']"}),
            'homologs': ('cyano.history.HistoryManyToManyField', [], {'blank': 'True', 'related_name': "u'genes'", 'null': 'True', 'symmetrical': 'False', 'to': u"orm['cyano.Homolog']"}),
            'is_essential': ('cyano.history.HistoryForeignKey', [], {'blank': 'True', 'related_name': "u'+'", 'null': 'True', 'to': u"orm['cyano.EntryBooleanData']"}),
            'length': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'parent_ptr_molecule': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "u'child_ptr_gene'", 'unique': 'True', 'primary_key': 'True', 'to': u"orm['cyano.Molecule']"}),
            'symbol': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255', 'blank': 'True'})
        },
        u'cyano.genome': {
            'Meta': {'object_name': 'Genome', '_ormbases': [u'cyano.Molecule']},
            'length': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'parent_ptr_molecule': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "u'child_ptr_genome'", 'unique': 'True', 'primary_key': 'True', 'to': u"orm['cyano.Molecule']"}),
            'sequence': ('django.db.models.fields.TextField', [], {'default': "u''", 'blank': 'True'})
        },
        u'cyano.grouppermission': {
            'Meta': {'unique_together': "((u'entry', u'group'),)", 'object_name': 'GroupPermission'},
            'allow': ('cyano.history.HistoryManyToManyField', [], {'related_name': "u'group_permission_allow'", 'symmetrical': 'False', 'to': u"orm['cyano.Permission']"}),
            'deny': ('cyano.history.HistoryManyToManyField', [], {'related_name': "u'group_permission_deny'", 'symmetrical': 'False', 'to': u"orm['cyano.Permission']"}),
            'entry': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'group_permissions'", 'to': u"orm['cyano.Entry']"}),
            'group': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'permissions'", 'to': u"orm['cyano.GroupProfile']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'})
        },
        u'cyano.groupprofile': {
            'Meta': {'object_name': 'GroupProfile'},
            'description': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255', 'blank': 'True'}),
            'group': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "u'profile'", 'unique': 'True', 'to': u"orm['auth.Group']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'})
        },
        u'cyano.homolog': {
            'Meta': {'ordering': "[u'xid']", 'object_name': 'Homolog'},
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'species': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'xid': ('django.db.models.fields.CharField', [], {'max_length': '255'})
        },
        u'cyano.kinetics': {
            'Meta': {'object_name': 'Kinetics'},
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'km': ('django.db.models.fields.CharField', [], {'max_length': '255', 'blank': 'True'}),
            'rate_law': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255', 'blank': 'True'}),
            'vmax': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'vmax_unit': ('django.db.models.fields.CharField', [], {'max_length': '255', 'blank': 'True'})
        },
        u'cyano.massspectrometryjob': {
            'Meta': {'object_name': 'MassSpectrometryJob', '_ormbases': [u'cyano.SpeciesComponent']},
            'parent_ptr_species_component': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "u'child_ptr_mass_spectrometry'", 'unique': 'True', 'primary_key': 'True', 'to': u"orm['cyano.SpeciesComponent']"})
        },
        u'cyano.massspectrometryprotein': {
            'Meta': {'object_name': 'MassSpectrometryProtein', '_ormbases': [u'cyano.Protein']},
            'ambiguous': ('cyano.history.HistoryManyToManyField', [], {'related_name': "u'ambiguous'", 'symmetrical': 'False', 'to': u"orm['cyano.EntryBasicTextData']"}),
            'coverage': ('django.db.models.fields.FloatField', [], {}),
            'length': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'mass': ('django.db.models.fields.FloatField', [], {}),
            'parent_ptr_species_component': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "u'child_ptr_ms_protein'", 'unique': 'True', 'to': u"orm['cyano.SpeciesComponent']"}),
            'pi': ('django.db.models.fields.FloatField', [], {}),
            u'protein_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': u"orm['cyano.Protein']", 'unique': 'True', 'primary_key': 'True'}),
            'score': ('django.db.models.fields.FloatField', [], {}),
            'sequence': ('django.db.models.fields.TextField', [], {}),
            'sub': ('cyano.history.HistoryManyToManyField', [], {'related_name': "u'sub'", 'symmetrical': 'False', 'to': u"orm['cyano.EntryBasicTextData']"})
        },
        u'cyano.massspectrometryproteindetail': {
            'Meta': {'object_name': 'MassSpectrometryProteinDetail'},
            'charge': ('django.db.models.fields.IntegerField', [], {}),
            'coordinate': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'delta_mass': ('django.db.models.fields.FloatField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'length': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'mass': ('django.db.models.fields.FloatField', [], {}),
            'missed_cleavages': ('django.db.models.fields.IntegerField', [], {}),
            'protein': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'protein_details'", 'to': u"orm['cyano.MassSpectrometryProtein']"}),
            'proteotypic': ('django.db.models.fields.BooleanField', [], {}),
            'retention_time': ('django.db.models.fields.IntegerField', [], {}),
            'sequence': ('django.db.models.fields.TextField', [], {}),
            'sequence_ptm': ('django.db.models.fields.TextField', [], {}),
            'theoretical_mass': ('django.db.models.fields.FloatField', [], {}),
            'zscore': ('django.db.models.fields.FloatField', [], {})
        },
        u'cyano.mediacomposition': {
            'Meta': {'ordering': "[u'-concentration']", 'object_name': 'MediaComposition'},
            'concentration': ('django.db.models.fields.FloatField', [], {}),
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'is_diffused': ('django.db.models.fields.BooleanField', [], {})
        },
        u'cyano.metabolite': {
            'Meta': {'object_name': 'Metabolite', '_ormbases': [u'cyano.Molecule']},
            'biomass_composition': ('cyano.history.HistoryManyToManyField', [], {'blank': 'True', 'related_name': "u'metabolites'", 'null': 'True', 'symmetrical': 'False', 'to': u"orm['cyano.BiomassComposition']"}),
            'charge': ('django.db.models.fields.IntegerField', [], {}),
            'deltag_formation': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'empirical_formula': ('django.db.models.fields.TextField', [], {}),
            'is_hydrophobic': ('django.db.models.fields.BooleanField', [], {}),
            'iupac_name': ('django.db.models.fields.TextField', [], {'default': "u''", 'blank': 'True'}),
            'log_d': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'log_p': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'map_coordinates': ('cyano.history.HistoryManyToManyField', [], {'blank': 'True', 'related_name': "u'metabolites'", 'null': 'True', 'symmetrical': 'False', 'to': u"orm['cyano.MetaboliteMapCoordinate']"}),
            'media_composition': ('cyano.history.HistoryForeignKey', [], {'blank': 'True', 'related_name': "u'metabolites'", 'null': 'True', 'on_delete': 'models.SET_NULL', 'to': u"orm['cyano.MediaComposition']"}),
            'parent_ptr_molecule': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "u'child_ptr_metabolite'", 'unique': 'True', 'primary_key': 'True', 'to': u"orm['cyano.Molecule']"}),
            'pi': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'pka': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'smiles': ('django.db.models.fields.TextField', [], {'default': "u''", 'blank': 'True'}),
            'traditional_name': ('django.db.models.fields.TextField', [], {'default': "u''", 'blank': 'True'}),
            'volume': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'})
        },
        u'cyano.metabolitemapcoordinate': {
            'Meta': {'ordering': "[u'x', u'y', u'compartment']", 'object_name': 'MetaboliteMapCoordinate'},
            'compartment': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.Compartment']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'x': ('django.db.models.fields.FloatField', [], {}),
            'y': ('django.db.models.fields.FloatField', [], {})
        },
        u'cyano.modificationreaction': {
            'Meta': {'object_name': 'ModificationReaction'},
            'compartment': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.Compartment']"}),
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'molecule': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'modification_reactions'", 'to': u"orm['cyano.Molecule']"}),
            'position': ('django.db.models.fields.PositiveIntegerField', [], {'null': 'True', 'blank': 'True'})
        },
        u'cyano.molecule': {
            'Meta': {'object_name': 'Molecule', '_ormbases': [u'cyano.SpeciesComponent']},
            'parent_ptr_species_component': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "u'child_ptr_molecule'", 'unique': 'True', 'primary_key': 'True', 'to': u"orm['cyano.SpeciesComponent']"})
        },
        u'cyano.note': {
            'Meta': {'object_name': 'Note', '_ormbases': [u'cyano.SpeciesComponent']},
            'parent_ptr_species_component': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "u'child_ptr_note'", 'unique': 'True', 'primary_key': 'True', 'to': u"orm['cyano.SpeciesComponent']"})
        },
        u'cyano.parameter': {
            'Meta': {'object_name': 'Parameter', '_ormbases': [u'cyano.SpeciesComponent']},
            'molecules': ('cyano.history.HistoryManyToManyField', [], {'blank': 'True', 'related_name': "u'parameters'", 'null': 'True', 'symmetrical': 'False', 'to': u"orm['cyano.Molecule']"}),
            'parent_ptr_species_component': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "u'child_ptr_parameter'", 'unique': 'True', 'primary_key': 'True', 'to': u"orm['cyano.SpeciesComponent']"}),
            'process': ('cyano.history.HistoryForeignKey', [], {'blank': 'True', 'related_name': "u'parameters'", 'null': 'True', 'on_delete': 'models.SET_NULL', 'to': u"orm['cyano.Process']"}),
            'reactions': ('cyano.history.HistoryManyToManyField', [], {'blank': 'True', 'related_name': "u'parameters'", 'null': 'True', 'symmetrical': 'False', 'to': u"orm['cyano.Reaction']"}),
            'state': ('cyano.history.HistoryForeignKey', [], {'blank': 'True', 'related_name': "u'parameters'", 'null': 'True', 'on_delete': 'models.SET_NULL', 'to': u"orm['cyano.State']"}),
            'value': ('cyano.history.HistoryForeignKey', [], {'to': u"orm['cyano.EntryCharData']"})
        },
        u'cyano.pathway': {
            'Meta': {'object_name': 'Pathway', '_ormbases': [u'cyano.SpeciesComponent']},
            'parent_ptr_species_component': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "u'child_ptr_pathway'", 'unique': 'True', 'primary_key': 'True', 'to': u"orm['cyano.SpeciesComponent']"})
        },
        u'cyano.peptide': {
            'Meta': {'object_name': 'Peptide', '_ormbases': [u'cyano.Protein']},
            'charge': ('django.db.models.fields.IntegerField', [], {}),
            'length': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'mass': ('django.db.models.fields.FloatField', [], {}),
            'parent_ptr_species_component': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "u'child_ptr_peptide'", 'unique': 'True', 'to': u"orm['cyano.SpeciesComponent']"}),
            u'protein_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': u"orm['cyano.Protein']", 'unique': 'True', 'primary_key': 'True'}),
            'proteins': ('cyano.history.HistoryManyToManyField', [], {'related_name': "u'peptides'", 'symmetrical': 'False', 'to': u"orm['cyano.EntryBasicTextData']"}),
            'proteotypic': ('django.db.models.fields.NullBooleanField', [], {'null': 'True', 'blank': 'True'}),
            'retention_time': ('django.db.models.fields.FloatField', [], {}),
            'sequence': ('django.db.models.fields.TextField', [], {'default': "u''", 'blank': 'True'}),
            'zscore': ('django.db.models.fields.FloatField', [], {})
        },
        u'cyano.permission': {
            'Meta': {'object_name': 'Permission'},
            'description': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255'})
        },
        u'cyano.plasmid': {
            'Meta': {'object_name': 'Plasmid', '_ormbases': [u'cyano.Genome']},
            'parent_ptr_genome': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "u'child_ptr_plasmid'", 'unique': 'True', 'primary_key': 'True', 'to': u"orm['cyano.Genome']"})
        },
        u'cyano.process': {
            'Meta': {'object_name': 'Process', '_ormbases': [u'cyano.SpeciesComponent']},
            'initialization_order': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'parent_ptr_species_component': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "u'child_ptr_process'", 'unique': 'True', 'primary_key': 'True', 'to': u"orm['cyano.SpeciesComponent']"})
        },
        u'cyano.prostheticgroupparticipant': {
            'Meta': {'object_name': 'ProstheticGroupParticipant'},
            'coefficient': ('django.db.models.fields.PositiveIntegerField', [], {'null': 'True', 'blank': 'True'}),
            'compartment': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.Compartment']"}),
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'metabolite': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'prosthetic_group_participants'", 'to': u"orm['cyano.Metabolite']"})
        },
        u'cyano.protein': {
            'Meta': {'object_name': 'Protein', '_ormbases': [u'cyano.Molecule']},
            'chaperones': ('cyano.history.HistoryManyToManyField', [], {'blank': 'True', 'related_name': "u'chaperone_substrates'", 'null': 'True', 'symmetrical': 'False', 'to': u"orm['cyano.Protein']"}),
            'dna_footprint': ('cyano.history.HistoryForeignKey', [], {'blank': 'True', 'related_name': "u'proteins'", 'null': 'True', 'to': u"orm['cyano.DNAFootprint']"}),
            'parent_ptr_molecule': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "u'child_ptr_protein'", 'unique': 'True', 'primary_key': 'True', 'to': u"orm['cyano.Molecule']"}),
            'prosthetic_groups': ('cyano.history.HistoryManyToManyField', [], {'blank': 'True', 'related_name': "u'proteins'", 'null': 'True', 'symmetrical': 'False', 'to': u"orm['cyano.ProstheticGroupParticipant']"}),
            'regulatory_rule': ('cyano.history.HistoryForeignKey', [], {'blank': 'True', 'related_name': "u'+'", 'null': 'True', 'on_delete': 'models.SET_NULL', 'to': u"orm['cyano.EntryCharData']"})
        },
        u'cyano.proteincomplex': {
            'Meta': {'object_name': 'ProteinComplex', '_ormbases': [u'cyano.Protein']},
            'biosynthesis': ('cyano.history.HistoryManyToManyField', [], {'related_name': "u'protein_complexes'", 'symmetrical': 'False', 'to': u"orm['cyano.ProteinComplexBiosythesisParticipant']"}),
            'disulfide_bonds': ('cyano.history.HistoryManyToManyField', [], {'blank': 'True', 'related_name': "u'protein_complexes'", 'null': 'True', 'symmetrical': 'False', 'to': u"orm['cyano.DisulfideBond']"}),
            'formation_process': ('cyano.history.HistoryForeignKey', [], {'blank': 'True', 'related_name': "u'formed_complexes'", 'null': 'True', 'on_delete': 'models.SET_NULL', 'to': u"orm['cyano.Process']"}),
            'parent_ptr_protein': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "u'child_ptr_protein_complex'", 'unique': 'True', 'primary_key': 'True', 'to': u"orm['cyano.Protein']"})
        },
        u'cyano.proteincomplexbiosythesisparticipant': {
            'Meta': {'object_name': 'ProteinComplexBiosythesisParticipant'},
            'coefficient': ('django.db.models.fields.FloatField', [], {}),
            'compartment': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.Compartment']"}),
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'molecule': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'protein_complex_biosythesis_participants'", 'to': u"orm['cyano.Molecule']"}),
            'residue': ('django.db.models.fields.PositiveIntegerField', [], {'null': 'True', 'blank': 'True'})
        },
        u'cyano.proteinmonomer': {
            'Meta': {'object_name': 'ProteinMonomer', '_ormbases': [u'cyano.Protein']},
            'gene': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'protein_monomers'", 'to': u"orm['cyano.Gene']"}),
            'is_n_terminal_methionine_cleaved': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'null': 'True', 'to': u"orm['cyano.EntryBooleanData']"}),
            'localization': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'protein_monomers'", 'null': 'True', 'to': u"orm['cyano.Compartment']"}),
            'parent_ptr_protein': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "u'child_ptr_protein_monomer'", 'unique': 'True', 'primary_key': 'True', 'to': u"orm['cyano.Protein']"}),
            'signal_sequence': ('cyano.history.HistoryForeignKey', [], {'blank': 'True', 'related_name': "u'protein_monomers'", 'null': 'True', 'on_delete': 'models.SET_NULL', 'to': u"orm['cyano.SignalSequence']"})
        },
        u'cyano.publicationreference': {
            'Meta': {'object_name': 'PublicationReference', '_ormbases': [u'cyano.SpeciesComponent']},
            'authors': ('django.db.models.fields.TextField', [], {'default': "u''", 'blank': 'True'}),
            'editors': ('django.db.models.fields.TextField', [], {'default': "u''", 'blank': 'True'}),
            'issue': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255', 'blank': 'True'}),
            'pages': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255', 'blank': 'True'}),
            'parent_ptr_species_component': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "u'child_ptr_publicationreference'", 'unique': 'True', 'primary_key': 'True', 'to': u"orm['cyano.SpeciesComponent']"}),
            'publication': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255', 'blank': 'True'}),
            'publisher': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255', 'blank': 'True'}),
            'title': ('django.db.models.fields.TextField', [], {'default': "u''", 'blank': 'True'}),
            'volume': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255', 'blank': 'True'}),
            'year': ('django.db.models.fields.PositiveIntegerField', [], {'null': 'True', 'blank': 'True'})
        },
        u'cyano.reaction': {
            'Meta': {'object_name': 'Reaction', '_ormbases': [u'cyano.SpeciesComponent']},
            'coenzymes': ('cyano.history.HistoryManyToManyField', [], {'blank': 'True', 'related_name': "u'reactions'", 'null': 'True', 'symmetrical': 'False', 'to': u"orm['cyano.CoenzymeParticipant']"}),
            'delta_g': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'direction': ('django.db.models.fields.CharField', [], {'max_length': '1'}),
            'enzyme': ('cyano.history.HistoryForeignKey', [], {'blank': 'True', 'related_name': "u'reactions'", 'null': 'True', 'on_delete': 'models.SET_NULL', 'to': u"orm['cyano.EnzymeParticipant']"}),
            'is_spontaneous': ('django.db.models.fields.BooleanField', [], {}),
            'keq': ('cyano.history.HistoryForeignKey', [], {'blank': 'True', 'related_name': "u'+'", 'null': 'True', 'on_delete': 'models.SET_NULL', 'to': u"orm['cyano.EntryPositiveFloatData']"}),
            'kinetics_backward': ('cyano.history.HistoryForeignKey', [], {'blank': 'True', 'related_name': "u'reactions_backward'", 'null': 'True', 'on_delete': 'models.SET_NULL', 'to': u"orm['cyano.Kinetics']"}),
            'kinetics_forward': ('cyano.history.HistoryForeignKey', [], {'blank': 'True', 'related_name': "u'reactions_forward'", 'null': 'True', 'on_delete': 'models.SET_NULL', 'to': u"orm['cyano.Kinetics']"}),
            'map_coordinates': ('cyano.history.HistoryManyToManyField', [], {'blank': 'True', 'related_name': "u'reactions'", 'null': 'True', 'symmetrical': 'False', 'to': u"orm['cyano.ReactionMapCoordinate']"}),
            'modification': ('cyano.history.HistoryForeignKey', [], {'blank': 'True', 'related_name': "u'reactions'", 'null': 'True', 'on_delete': 'models.SET_NULL', 'to': u"orm['cyano.ModificationReaction']"}),
            'optimal_ph': ('cyano.history.HistoryForeignKey', [], {'blank': 'True', 'related_name': "u'+'", 'null': 'True', 'on_delete': 'models.SET_NULL', 'to': u"orm['cyano.EntryPositiveFloatData']"}),
            'optimal_temperature': ('cyano.history.HistoryForeignKey', [], {'blank': 'True', 'related_name': "u'+'", 'null': 'True', 'on_delete': 'models.SET_NULL', 'to': u"orm['cyano.EntryFloatData']"}),
            'parent_ptr_species_component': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "u'child_ptr_reaction'", 'unique': 'True', 'primary_key': 'True', 'to': u"orm['cyano.SpeciesComponent']"}),
            'pathways': ('cyano.history.HistoryManyToManyField', [], {'blank': 'True', 'related_name': "u'reactions'", 'null': 'True', 'symmetrical': 'False', 'to': u"orm['cyano.Pathway']"}),
            'processes': ('cyano.history.HistoryForeignKey', [], {'blank': 'True', 'related_name': "u'reactions'", 'null': 'True', 'on_delete': 'models.SET_NULL', 'to': u"orm['cyano.Process']"}),
            'states': ('cyano.history.HistoryForeignKey', [], {'blank': 'True', 'related_name': "u'reactions'", 'null': 'True', 'on_delete': 'models.SET_NULL', 'to': u"orm['cyano.State']"}),
            'stoichiometry': ('cyano.history.HistoryManyToManyField', [], {'related_name': "u'reactions'", 'symmetrical': 'False', 'to': u"orm['cyano.ReactionStoichiometryParticipant']"})
        },
        u'cyano.reactionmapcoordinate': {
            'Meta': {'ordering': "[u'value_x', u'value_y', u'label_x', u'label_y']", 'object_name': 'ReactionMapCoordinate'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'label_x': ('django.db.models.fields.FloatField', [], {}),
            'label_y': ('django.db.models.fields.FloatField', [], {}),
            'path': ('django.db.models.fields.TextField', [], {}),
            'value_x': ('django.db.models.fields.FloatField', [], {}),
            'value_y': ('django.db.models.fields.FloatField', [], {})
        },
        u'cyano.reactionstoichiometryparticipant': {
            'Meta': {'object_name': 'ReactionStoichiometryParticipant'},
            'coefficient': ('django.db.models.fields.FloatField', [], {}),
            'compartment': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.Compartment']"}),
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'molecule': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'reaction_stoichiometry_participants'", 'to': u"orm['cyano.Molecule']"})
        },
        u'cyano.revision': {
            'Meta': {'object_name': 'Revision'},
            'action': ('django.db.models.fields.CharField', [], {'max_length': '1'}),
            'content_type': ('cyano.history.HistoryForeignKey', [], {'to': u"orm['contenttypes.ContentType']"}),
            'detail': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'revisions'", 'to': u"orm['cyano.RevisionDetail']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'new_data': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'object_id': ('django.db.models.fields.IntegerField', [], {'db_index': 'True'})
        },
        u'cyano.revisiondetail': {
            'Meta': {'object_name': 'RevisionDetail'},
            'date': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'reason': ('django.db.models.fields.TextField', [], {'default': "u''", 'blank': 'True'}),
            'user': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.UserProfile']"})
        },
        u'cyano.signalsequence': {
            'Meta': {'ordering': "[u'type', u'location', u'length']", 'object_name': 'SignalSequence'},
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'length': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'location': ('django.db.models.fields.CharField', [], {'max_length': '1'}),
            'type': ('django.db.models.fields.CharField', [], {'max_length': '20'})
        },
        u'cyano.species': {
            'Meta': {'ordering': "[u'wid']", 'object_name': 'Species', '_ormbases': [u'cyano.Entry']},
            'cross_references': ('cyano.history.HistoryManyToManyField', [], {'blank': 'True', 'related_name': "u'cross_referenced_species'", 'null': 'True', 'symmetrical': 'False', 'to': u"orm['cyano.CrossReference']"}),
            'genetic_code': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'parent_ptr_entry': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "u'child_ptr_species'", 'unique': 'True', 'primary_key': 'True', 'to': u"orm['cyano.Entry']"}),
            'publication_references': ('cyano.history.HistoryManyToManyField', [], {'blank': 'True', 'related_name': "u'publication_referenced_species'", 'null': 'True', 'symmetrical': 'False', 'to': u"orm['cyano.PublicationReference']"})
        },
        u'cyano.speciescomponent': {
            'Meta': {'object_name': 'SpeciesComponent'},
            'cross_references': ('cyano.history.HistoryManyToManyField', [], {'blank': 'True', 'related_name': "u'cross_referenced_components'", 'null': 'True', 'symmetrical': 'False', 'to': u"orm['cyano.CrossReference']"}),
            u'entry_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': u"orm['cyano.Entry']", 'unique': 'True', 'primary_key': 'True'}),
            'parent': ('cyano.history.HistoryForeignKey', [], {'blank': 'True', 'related_name': "u'children'", 'null': 'True', 'on_delete': 'models.SET_NULL', 'to': u"orm['cyano.SpeciesComponent']"}),
            'publication_references': ('cyano.history.HistoryManyToManyField', [], {'blank': 'True', 'related_name': "u'publication_referenced_components'", 'null': 'True', 'symmetrical': 'False', 'to': u"orm['cyano.PublicationReference']"}),
            'species': ('cyano.history.HistoryManyToManyField', [], {'related_name': "u'species_components'", 'symmetrical': 'False', 'to': u"orm['cyano.Species']"}),
            'type': ('cyano.history.HistoryManyToManyField', [], {'blank': 'True', 'related_name': "u'members'", 'null': 'True', 'symmetrical': 'False', 'to': u"orm['cyano.Type']"})
        },
        u'cyano.state': {
            'Meta': {'object_name': 'State', '_ormbases': [u'cyano.SpeciesComponent']},
            'parent_ptr_species_component': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "u'child_ptr_state'", 'unique': 'True', 'primary_key': 'True', 'to': u"orm['cyano.SpeciesComponent']"})
        },
        u'cyano.stimulus': {
            'Meta': {'object_name': 'Stimulus', '_ormbases': [u'cyano.Molecule']},
            'parent_ptr_molecule': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "u'child_ptr_stimulus'", 'unique': 'True', 'primary_key': 'True', 'to': u"orm['cyano.Molecule']"}),
            'value': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.EntryFloatData']"})
        },
        u'cyano.synonym': {
            'Meta': {'ordering': "[u'name']", 'object_name': 'Synonym'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '255'})
        },
        u'cyano.tablemeta': {
            'Meta': {'object_name': 'TableMeta'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'model_name': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '255'}),
            'table_name': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '255'})
        },
        u'cyano.transcriptionalregulation': {
            'Meta': {'object_name': 'TranscriptionalRegulation', '_ormbases': [u'cyano.SpeciesComponent']},
            'activity': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.EntryPositiveFloatData']"}),
            'affinity': ('cyano.history.HistoryForeignKey', [], {'blank': 'True', 'related_name': "u'+'", 'null': 'True', 'on_delete': 'models.SET_NULL', 'to': u"orm['cyano.EntryPositiveFloatData']"}),
            'binding_site': ('cyano.history.HistoryForeignKey', [], {'blank': 'True', 'related_name': "u'transcriptional_regulations'", 'null': 'True', 'on_delete': 'models.SET_NULL', 'to': u"orm['cyano.BindingSite']"}),
            'parent_ptr_species_component': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "u'child_ptr_transcriptional_regulation'", 'unique': 'True', 'primary_key': 'True', 'to': u"orm['cyano.SpeciesComponent']"}),
            'transcription_factor': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'transcriptional_regulations'", 'to': u"orm['cyano.Protein']"}),
            'transcription_unit': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'transcriptional_regulations'", 'to': u"orm['cyano.TranscriptionUnit']"})
        },
        u'cyano.transcriptionunit': {
            'Meta': {'object_name': 'TranscriptionUnit', '_ormbases': [u'cyano.Molecule']},
            'genes': ('cyano.history.HistoryManyToManyField', [], {'related_name': "u'transcription_units'", 'symmetrical': 'False', 'to': u"orm['cyano.Gene']"}),
            'parent_ptr_molecule': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "u'child_ptr_transcription_unit'", 'unique': 'True', 'primary_key': 'True', 'to': u"orm['cyano.Molecule']"}),
            'promoter_10_coordinate': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'promoter_10_length': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'promoter_35_coordinate': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'promoter_35_length': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'tss_coordinate': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'})
        },
        u'cyano.type': {
            'Meta': {'object_name': 'Type', '_ormbases': [u'cyano.SpeciesComponent']},
            'parent_ptr_species_component': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "u'child_ptr_type'", 'unique': 'True', 'primary_key': 'True', 'to': u"orm['cyano.SpeciesComponent']"})
        },
        u'cyano.userpermission': {
            'Meta': {'unique_together': "((u'entry', u'user'),)", 'object_name': 'UserPermission'},
            'allow': ('cyano.history.HistoryManyToManyField', [], {'related_name': "u'user_permission_allow'", 'symmetrical': 'False', 'to': u"orm['cyano.Permission']"}),
            'deny': ('cyano.history.HistoryManyToManyField', [], {'related_name': "u'user_permission_deny'", 'symmetrical': 'False', 'to': u"orm['cyano.Permission']"}),
            'entry': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'user_permissions'", 'to': u"orm['cyano.Entry']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'user': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'permissions'", 'to': u"orm['cyano.UserProfile']"})
        },
        u'cyano.userprofile': {
            'Meta': {'ordering': "[u'user__last_name', u'user__first_name']", 'object_name': 'UserProfile'},
            'address': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255', 'blank': 'True'}),
            'affiliation': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255', 'blank': 'True'}),
            'city': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255', 'blank': 'True'}),
            'country': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255', 'blank': 'True'}),
            'force_password_change': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'phone': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255', 'blank': 'True'}),
            'state': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255', 'blank': 'True'}),
            'user': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "u'profile'", 'unique': 'True', 'to': u"orm['auth.User']"}),
            'website': ('django.db.models.fields.URLField', [], {'default': "u''", 'max_length': '255', 'blank': 'True'}),
            'zip': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255', 'blank': 'True'})
        }
    }

    complete_apps = ['cyano']