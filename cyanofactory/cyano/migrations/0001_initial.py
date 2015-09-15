# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import datetime
import re
import cyano.models
import cyano.history
import django.db.models.deletion
from django.conf import settings
import django.core.validators


class Migration(migrations.Migration):

    dependencies = [
        ('contenttypes', '0002_remove_content_type_name'),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('auth', '0006_require_contenttypes_0002'),
    ]

    operations = [
        migrations.CreateModel(
            name='Basket',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(default='', max_length=255, verbose_name='Basket name')),
            ],
        ),
        migrations.CreateModel(
            name='BasketComponent',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('basket', cyano.history.HistoryForeignKey(related_name='components', verbose_name='In basket', to='cyano.Basket')),
            ],
        ),
        migrations.CreateModel(
            name='BindingSite',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('coordinate', models.PositiveIntegerField(verbose_name='Coordinate (nt)')),
                ('length', models.PositiveIntegerField(verbose_name='Length (nt)')),
                ('direction', models.CharField(max_length=10, verbose_name='Direction', choices=[('f', 'Forward'), ('r', 'Reverse')])),
            ],
            options={
                'ordering': ['coordinate', 'length'],
                'verbose_name': 'Binding site',
                'verbose_name_plural': 'Binding sites',
            },
        ),
        migrations.CreateModel(
            name='BiomassComposition',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('concentration', models.FloatField(verbose_name='Concentration (mmol gDCW<sup>-1</sup>)', validators=[django.core.validators.MinValueValidator(0)])),
            ],
            options={
                'ordering': ['-concentration'],
                'verbose_name': 'Biomass composition',
                'verbose_name_plural': 'Biomass composition',
            },
        ),
        migrations.CreateModel(
            name='Codon',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('sequence', models.CharField(max_length=3, verbose_name='Sequence', validators=[cyano.models.validate_dna_sequence, django.core.validators.MinLengthValidator(3)])),
            ],
            options={
                'ordering': ['sequence'],
                'verbose_name': 'Codon',
                'verbose_name_plural': 'Codons',
            },
        ),
        migrations.CreateModel(
            name='CoenzymeParticipant',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('coefficient', models.FloatField(blank=True, null=True, verbose_name='Coefficient', validators=[django.core.validators.MinValueValidator(0)])),
            ],
            options={
                'ordering': [],
                'verbose_name': 'Coenzyme participant',
                'verbose_name_plural': 'Coenzyme participants',
            },
        ),
        migrations.CreateModel(
            name='CrossReference',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('xid', models.CharField(max_length=255, verbose_name='External ID')),
                ('source', models.CharField(max_length=20, verbose_name='Source', choices=[('ATCC', 'ATCC'), ('BiGG', 'BiGG'), ('BioCyc', 'BioCyc'), ('BioProject', 'BioProject'), ('CAS', 'CAS'), ('ChEBI', 'ChEBI'), ('CMR', 'CMR'), ('EC', 'EC'), ('GenBank', 'GenBank'), ('ISBN', 'ISBN'), ('KEGG', 'KEGG'), ('KNApSAcK', 'KNApSAcK'), ('LipidBank', 'LipidBank'), ('LIPIDMAPS', 'LIPIDMAPS'), ('PDB', 'PDB'), ('PDBCCD', 'PDBCCD'), ('PubChem', 'PubChem'), ('PubMed', 'PubMed'), ('RefSeq', 'RefSeq'), ('SABIO-RK', 'SABIO-RK'), ('SwissProt', 'SwissProt'), ('Taxonomy', 'Taxonomy'), ('ThreeDMET', 'ThreeDMET'), ('URL', 'URL')])),
            ],
            options={
                'ordering': ['xid'],
                'verbose_name': 'Cross reference',
                'verbose_name_plural': 'Cross references',
            },
        ),
        migrations.CreateModel(
            name='DisulfideBond',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('residue_1', models.PositiveIntegerField(verbose_name='Residue-1')),
                ('residue_2', models.PositiveIntegerField(verbose_name='Residue-2')),
            ],
            options={
                'ordering': [],
                'verbose_name': 'Disulfide bond',
                'verbose_name_plural': 'Disulfide bonds',
            },
        ),
        migrations.CreateModel(
            name='DNAFootprint',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('length', models.PositiveIntegerField(null=True, verbose_name='Length (nt)', blank=True)),
                ('binding', models.CharField(default='', max_length=10, verbose_name='Binding', blank=True, choices=[('dsDNA', 'dsDNA'), ('ssDNA', 'ssDNA'), ('xsDNA', 'xsDNA')])),
                ('region', models.CharField(default='', max_length=10, verbose_name='Region', blank=True, choices=[('dsDNA', 'dsDNA'), ('ssDNA', 'ssDNA'), ('xsDNA', 'xsDNA')])),
            ],
            options={
                'ordering': [],
                'verbose_name': 'DNA footprint',
                'verbose_name_plural': 'DNA footprints',
            },
        ),
        migrations.CreateModel(
            name='Entry',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('wid', models.SlugField(max_length=150, verbose_name='WID', validators=[django.core.validators.RegexValidator(re.compile('^[-a-zA-Z0-9_]+\\Z'), "Enter a valid 'slug' consisting of letters, numbers, underscores or hyphens.", 'invalid')])),
                ('name', models.CharField(default='', max_length=255, verbose_name='Name', blank=True)),
                ('comments', models.TextField(default='', verbose_name='Comments', blank=True)),
            ],
            options={
                'wid_unique': False,
                'concrete_entry_model': False,
                'get_latest_by': 'createdDate',
                'ordering': ['wid'],
                'verbose_name_plural': 'Entries',
                'facet_fields': [],
                'listing': ['wid', 'name'],
                'verbose_name': 'Entry',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'comments'],
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Comments', {'fields': ['comments']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
                'permissions': (('view_normal', 'View entry'), ('view_delete', 'View deleted revisions of entry'), ('view_permission', 'View permissions of entry'), ('view_history', 'View older (not deleted) revisions of entry'), ('edit_permission', 'Allow modifying of permissions')),
            },
        ),
        migrations.CreateModel(
            name='EntryBasicTextData',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('value', models.TextField(default='', verbose_name='Value')),
            ],
            options={
                'verbose_name': 'Entry basic text data',
                'verbose_name_plural': 'Entry basic text data',
            },
        ),
        migrations.CreateModel(
            name='EntryBooleanData',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('value', models.BooleanField(verbose_name='Value')),
            ],
            options={
                'ordering': ['value'],
                'verbose_name': 'Entry Boolean data',
                'verbose_name_plural': 'Entry Boolean data',
            },
        ),
        migrations.CreateModel(
            name='EntryCharData',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('value', models.CharField(default='', max_length=255, verbose_name='Value', blank=True)),
                ('units', models.CharField(default='', max_length=255, verbose_name='Units', blank=True)),
            ],
            options={
                'ordering': ['value', 'units'],
                'verbose_name': 'Entry char data',
                'verbose_name_plural': 'Entry char data',
            },
        ),
        migrations.CreateModel(
            name='EntryFloatData',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('value', models.FloatField(verbose_name='Value')),
                ('units', models.CharField(default='', max_length=255, verbose_name='Units', blank=True)),
            ],
            options={
                'ordering': ['value', 'units'],
                'verbose_name': 'Entry float data',
                'verbose_name_plural': 'Entry float data',
            },
        ),
        migrations.CreateModel(
            name='EntryGroupObjectPermission',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='EntryPositiveFloatData',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('value', models.FloatField(verbose_name='Value', validators=[django.core.validators.MinValueValidator(0)])),
                ('units', models.CharField(default='', max_length=255, verbose_name='Units', blank=True)),
            ],
            options={
                'ordering': ['value', 'units'],
                'verbose_name': 'Entry positive float data',
                'verbose_name_plural': 'Entry positive float data',
            },
        ),
        migrations.CreateModel(
            name='EntryTextData',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('value', models.TextField(default='', verbose_name='Value', blank=True)),
                ('units', models.CharField(default='', max_length=255, verbose_name='Units', blank=True)),
            ],
            options={
                'ordering': ['value', 'units'],
                'verbose_name': 'Entry text data',
                'verbose_name_plural': 'Entry text data',
            },
        ),
        migrations.CreateModel(
            name='EntryUserObjectPermission',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='EnzymeParticipant',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
            ],
            options={
                'ordering': [],
                'verbose_name': 'Enzyme participant',
                'verbose_name_plural': 'Enzyme participants',
            },
        ),
        migrations.CreateModel(
            name='Evidence',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('value', models.TextField(default='', verbose_name='Value', blank=True)),
                ('units', models.CharField(max_length=255, null=True, verbose_name='Units', blank=True)),
                ('is_experimentally_constrained', models.BooleanField(verbose_name='Is experimentally <br/>constrained')),
                ('species', models.CharField(max_length=255, null=True, verbose_name='Species', blank=True)),
                ('media', models.CharField(max_length=255, null=True, verbose_name='Media', blank=True)),
                ('pH', models.FloatField(null=True, verbose_name='pH', blank=True)),
                ('temperature', models.FloatField(null=True, verbose_name='Temperature (C)', blank=True)),
                ('comments', models.TextField(default='', verbose_name='Comments', blank=True)),
            ],
            options={
                'ordering': ['value', 'units'],
                'verbose_name': 'Evidence',
                'verbose_name_plural': 'Evidence',
            },
        ),
        migrations.CreateModel(
            name='FeaturePosition',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('coordinate', models.PositiveIntegerField(verbose_name='Coordinate (nt)')),
                ('length', models.PositiveIntegerField(verbose_name='Length (nt)')),
                ('direction', models.CharField(max_length=10, verbose_name='Direction', choices=[('f', 'Forward'), ('r', 'Reverse')])),
            ],
            options={
                'verbose_name': 'Feature Position',
                'verbose_name_plural': 'Feature Positions',
                'fieldsets': [('Structure', {'fields': [{'verbose_name': 'Structure', 'name': 'structure'}, {'verbose_name': 'Structure Filter', 'name': 'structure_filter'}, {'verbose_name': 'Sequence', 'name': 'sequence'}, {'verbose_name': 'Genes', 'name': 'genes'}, {'verbose_name': 'Transcription units', 'name': 'transcription_units'}]})],
                'field_list': ['id', 'chromosome_feature', 'chromosome', 'coordinate', 'length'],
            },
        ),
        migrations.CreateModel(
            name='GlobalPermission',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
            ],
            options={
                'permissions': (('access_species', 'Can access any species'), ('create_mutant', 'Can create mutants'), ('access_sbgn', 'Can access SBGN map')),
            },
        ),
        migrations.CreateModel(
            name='GroupProfile',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('description', models.CharField(default='', max_length=255, verbose_name='Group description', blank=True)),
                ('group', models.OneToOneField(related_name='profile', to='auth.Group')),
            ],
        ),
        migrations.CreateModel(
            name='Homolog',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('xid', models.CharField(max_length=255, verbose_name='External ID')),
                ('species', models.CharField(max_length=20, verbose_name='Species', choices=[('B. subtilis', 'B. subtilis'), ('E. coli', 'E. coli'), ('M. hyopneumoniae', 'M. hyopneumoniae'), ('M. mobile', 'M. mobile'), ('M. pneumoniae', 'M. pneumoniae'), ('S. coelicolor', 'S. coelicolor'), ('S. oneidensis', 'S. oneidensis')])),
                ('evidence', cyano.history.HistoryManyToManyField(to='cyano.Evidence', verbose_name='Evidence', blank=True)),
            ],
            options={
                'ordering': ['xid'],
                'verbose_name': 'Homolog',
                'verbose_name_plural': 'Homologs',
            },
        ),
        migrations.CreateModel(
            name='Kinetics',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('rate_law', models.CharField(default='', max_length=255, verbose_name='Rate law', blank=True)),
                ('km', models.CharField(blank=True, max_length=255, verbose_name='K<sub>m</sub> (&mu;M)', validators=[django.core.validators.RegexValidator('^([0-9\\.]+)(, [0-9\\.]+)*$')])),
                ('vmax', models.FloatField(blank=True, null=True, verbose_name='V<sub>max</sub>', validators=[django.core.validators.MinValueValidator(0)])),
                ('vmax_unit', models.CharField(blank=True, max_length=255, verbose_name='V<sub>max</sub> Unit', choices=[('1/min', '1/min'), ('U/mg', 'U/mg')])),
                ('evidence', cyano.history.HistoryManyToManyField(to='cyano.Evidence', verbose_name='Evidence', blank=True)),
            ],
            options={
                'ordering': [],
                'verbose_name': 'Kinetics',
                'verbose_name_plural': 'Kinetics',
            },
        ),
        migrations.CreateModel(
            name='MassSpectrometryProteinDetail',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('sequence', models.TextField(verbose_name='Sequence')),
                ('sequence_ptm', models.TextField(verbose_name='Sequence + PTMs')),
                ('coordinate', models.PositiveIntegerField(verbose_name='Coordinate (nt)')),
                ('length', models.PositiveIntegerField(verbose_name='Length (nt)')),
                ('proteotypic', models.BooleanField(verbose_name='Proteotypic')),
                ('zscore', models.FloatField(verbose_name='zscore')),
                ('delta_mass', models.FloatField(verbose_name='Delta Mass (ppm)')),
                ('mass', models.FloatField(verbose_name='Experimental Mass (m/z)')),
                ('charge', models.IntegerField(verbose_name='Charge')),
                ('retention_time', models.IntegerField(verbose_name='Retention Time (min)')),
                ('theoretical_mass', models.FloatField(verbose_name='Theoretical Mass (Da)')),
                ('missed_cleavages', models.IntegerField(verbose_name='Missed Cleavages')),
            ],
            options={
                'verbose_name': 'Mass Spectrometry Protein Detail',
                'verbose_name_plural': 'Mass Spectrometry Protein Details',
                'fieldsets': [('Mass Spectrometry Details', {'fields': ['sequence', 'sequence_ptm', 'proteotypic', 'zscore', 'delta_mass', 'mass', 'charge', 'retention_time', 'theoretical_mass', 'missed_cleavages']})],
                'field_list': ['id', 'chromosome_feature', 'chromosome', 'coordinate', 'length'],
            },
        ),
        migrations.CreateModel(
            name='MediaComposition',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('concentration', models.FloatField(verbose_name='Concentration (mM)', validators=[django.core.validators.MinValueValidator(0)])),
                ('is_diffused', models.BooleanField(verbose_name='Is diffused')),
                ('evidence', cyano.history.HistoryManyToManyField(to='cyano.Evidence', verbose_name='Evidence', blank=True)),
            ],
            options={
                'ordering': ['-concentration'],
                'verbose_name': 'Media composition',
                'verbose_name_plural': 'Media composition',
            },
        ),
        migrations.CreateModel(
            name='MetaboliteMapCoordinate',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('x', models.FloatField(verbose_name='X')),
                ('y', models.FloatField(verbose_name='Y')),
            ],
            options={
                'ordering': ['x', 'y', 'compartment'],
                'verbose_name': 'Metabolite map coordinate',
                'verbose_name_plural': 'Metabolite map coordinates',
            },
        ),
        migrations.CreateModel(
            name='ModificationReaction',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('position', models.PositiveIntegerField(null=True, verbose_name='Position', blank=True)),
                ('evidence', cyano.history.HistoryManyToManyField(to='cyano.Evidence', verbose_name='Evidence', blank=True)),
            ],
            options={
                'ordering': [],
                'verbose_name': 'Protein monomer',
                'verbose_name_plural': 'Protein monomers',
            },
        ),
        migrations.CreateModel(
            name='ProstheticGroupParticipant',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('coefficient', models.PositiveIntegerField(null=True, verbose_name='Coefficient', blank=True)),
                ('evidence', cyano.history.HistoryManyToManyField(to='cyano.Evidence', verbose_name='Evidence', blank=True)),
            ],
            options={
                'ordering': [],
                'verbose_name': 'Prosthetic group participant',
                'verbose_name_plural': 'Prosthetic group participants',
            },
        ),
        migrations.CreateModel(
            name='ProteinComparison',
            fields=[
                ('wid', models.AutoField(serialize=False, primary_key=True)),
                ('protein_a', models.IntegerField()),
                ('protein_b', models.IntegerField()),
                ('equal_score', models.FloatField()),
            ],
        ),
        migrations.CreateModel(
            name='ProteinComplexBiosythesisParticipant',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('residue', models.PositiveIntegerField(null=True, verbose_name='Residue', blank=True)),
                ('coefficient', models.FloatField(verbose_name='Coefficient')),
                ('evidence', cyano.history.HistoryManyToManyField(to='cyano.Evidence', verbose_name='Evidence', blank=True)),
            ],
            options={
                'ordering': [],
                'verbose_name': 'Protein complex biosythesis participant',
                'verbose_name_plural': 'Protein complex biosythesis participants',
            },
        ),
        migrations.CreateModel(
            name='ReactionMapCoordinate',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('path', models.TextField(verbose_name='Path')),
                ('value_x', models.FloatField(verbose_name='Label-X')),
                ('value_y', models.FloatField(verbose_name='Label-Y')),
                ('label_x', models.FloatField(verbose_name='Value-X')),
                ('label_y', models.FloatField(verbose_name='Value-Y')),
            ],
            options={
                'ordering': ['value_x', 'value_y', 'label_x', 'label_y'],
                'verbose_name': 'Reaction map coordinate',
                'verbose_name_plural': 'Reaction map coordinates',
            },
        ),
        migrations.CreateModel(
            name='ReactionStoichiometryParticipant',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('coefficient', models.FloatField(verbose_name='Coefficient')),
                ('evidence', cyano.history.HistoryManyToManyField(to='cyano.Evidence', verbose_name='Evidence', blank=True)),
            ],
            options={
                'ordering': [],
                'verbose_name': 'Molecule coefficient compartment',
                'verbose_name_plural': 'Molecule coefficient compartments',
            },
        ),
        migrations.CreateModel(
            name='Revision',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('object_id', models.IntegerField(verbose_name='Current version primary key', db_index=True)),
                ('action', models.CharField(max_length=1, choices=[('I', 'Insert'), ('U', 'Update'), ('D', 'Delete'), ('X', 'Unknown')])),
                ('new_data', models.TextField(blank=True)),
                ('content_type', cyano.history.HistoryForeignKey(to='contenttypes.ContentType')),
            ],
        ),
        migrations.CreateModel(
            name='RevisionDetail',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('date', models.DateTimeField(default=datetime.datetime.now, verbose_name='Modification date')),
                ('reason', models.TextField(default='', verbose_name='Reason for edit', blank=True)),
            ],
        ),
        migrations.CreateModel(
            name='SignalSequence',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('type', models.CharField(max_length=20, verbose_name='Type', choices=[('lipoprotein', 'Lipoprotein'), ('secretory', 'Secretory')])),
                ('location', models.CharField(max_length=1, verbose_name='Location', choices=[('N', 'N-terminus'), ('C', 'C-terminus')])),
                ('length', models.PositiveIntegerField(verbose_name='Length (nt)')),
                ('evidence', cyano.history.HistoryManyToManyField(to='cyano.Evidence', verbose_name='Evidence', blank=True)),
            ],
            options={
                'ordering': ['type', 'location', 'length'],
                'verbose_name': 'Signal sequence',
                'verbose_name_plural': 'Signal sequences',
            },
        ),
        migrations.CreateModel(
            name='Synonym',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=255, verbose_name='Name')),
            ],
            options={
                'ordering': ['name'],
                'verbose_name': 'Synonym',
                'verbose_name_plural': 'Synonyms',
            },
        ),
        migrations.CreateModel(
            name='TableMeta',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('table_name', models.CharField(unique=True, max_length=255, verbose_name='Name of the table')),
                ('model_name', models.CharField(unique=True, max_length=255, verbose_name='Name of the model associated with the table')),
            ],
        ),
        migrations.CreateModel(
            name='UserProfile',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('affiliation', models.CharField(default='', max_length=255, verbose_name='Affiliation', blank=True)),
                ('website', models.URLField(default='', max_length=255, verbose_name='Website', blank=True)),
                ('phone', models.CharField(default='', max_length=255, verbose_name='Phone', blank=True)),
                ('address', models.CharField(default='', max_length=255, verbose_name='Address', blank=True)),
                ('city', models.CharField(default='', max_length=255, verbose_name='City', blank=True)),
                ('state', models.CharField(default='', max_length=255, verbose_name='State', blank=True)),
                ('zip', models.CharField(default='', max_length=255, verbose_name='Zip', blank=True)),
                ('country', models.CharField(default='', max_length=255, verbose_name='Country', blank=True)),
                ('force_password_change', models.BooleanField(default=False)),
                ('user', models.OneToOneField(related_name='profile', to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'ordering': ['user__last_name', 'user__first_name'],
                'get_latest_by': 'user__date_joined',
                'verbose_name': 'User profile',
                'verbose_name_plural': 'User profiles',
            },
        ),
        migrations.CreateModel(
            name='Species',
            fields=[
                ('parent_ptr_entry', models.OneToOneField(parent_link=True, related_name='child_ptr_species', primary_key=True, serialize=False, to='cyano.Entry', verbose_name='Entry')),
                ('genetic_code', models.CharField(max_length=50, verbose_name='Genetic code', choices=[('1', 'Standard'), ('2', 'Vertebrate'), ('3', 'Yeast'), ('4', 'Mold, protozoa, coelenterate mitochondria, mycoplasma, and spiroplasma'), ('5', 'Invertebrate mitochondria'), ('6', 'Ciliate, dasycladacean and hexamita'), ('9', 'Echinoderm and flatworm mitochondria'), ('10', 'Euplotid'), ('11', 'Bacteria, archaea and plant plastids'), ('12', 'Alternative yeast'), ('13', 'Ascidian mitochondria'), ('14', 'Alternative flatworm mitochondria'), ('15', 'Blepharisma'), ('16', 'Chlorophycean mitochondria'), ('21', 'Trematode mitochondria'), ('22', 'Scenedesmus obliquus mitochondria'), ('23', 'Thraustochytrium mitochondria'), ('24', 'Pterobranchia mitochondria')])),
            ],
            options={
                'wid_unique': True,
                'concrete_entry_model': True,
                'facet_fields': [],
                'verbose_name_plural': 'Species',
                'verbose_name': 'Species',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'publication_references', 'cross_references', 'genetic_code', 'comments'],
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Physiology', {'fields': ['genetic_code']}), ('Comments', {'fields': ['comments']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
            },
            bases=('cyano.entry',),
        ),
        migrations.CreateModel(
            name='SpeciesComponent',
            fields=[
                ('entry_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='cyano.Entry')),
            ],
            options={
                'wid_unique': False,
                'concrete_entry_model': False,
                'facet_fields': ['type'],
                'verbose_name_plural': 'Species components',
                'verbose_name': 'Species component',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'comments', 'publication_references'],
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
            },
            bases=('cyano.entry',),
        ),
        migrations.AddField(
            model_name='revisiondetail',
            name='user',
            field=cyano.history.HistoryForeignKey(related_name='+', editable=False, to='cyano.UserProfile', verbose_name='Modified by'),
        ),
        migrations.AddField(
            model_name='revision',
            name='detail',
            field=cyano.history.HistoryForeignKey(related_name='revisions', editable=False, to='cyano.RevisionDetail', verbose_name='Details about this revision'),
        ),
        migrations.AddField(
            model_name='enzymeparticipant',
            name='evidence',
            field=cyano.history.HistoryManyToManyField(to='cyano.Evidence', verbose_name='Evidence', blank=True),
        ),
        migrations.AddField(
            model_name='entryuserobjectpermission',
            name='content_object',
            field=models.ForeignKey(to='cyano.Entry'),
        ),
        migrations.AddField(
            model_name='entryuserobjectpermission',
            name='permission',
            field=models.ForeignKey(to='auth.Permission'),
        ),
        migrations.AddField(
            model_name='entryuserobjectpermission',
            name='user',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL),
        ),
        migrations.AddField(
            model_name='entrytextdata',
            name='evidence',
            field=cyano.history.HistoryManyToManyField(to='cyano.Evidence', verbose_name='Evidence', blank=True),
        ),
        migrations.AddField(
            model_name='entrypositivefloatdata',
            name='evidence',
            field=cyano.history.HistoryManyToManyField(to='cyano.Evidence', verbose_name='Evidence', blank=True),
        ),
        migrations.AddField(
            model_name='entrygroupobjectpermission',
            name='content_object',
            field=models.ForeignKey(to='cyano.Entry'),
        ),
        migrations.AddField(
            model_name='entrygroupobjectpermission',
            name='group',
            field=models.ForeignKey(to='auth.Group'),
        ),
        migrations.AddField(
            model_name='entrygroupobjectpermission',
            name='permission',
            field=models.ForeignKey(to='auth.Permission'),
        ),
        migrations.AddField(
            model_name='entryfloatdata',
            name='evidence',
            field=cyano.history.HistoryManyToManyField(to='cyano.Evidence', verbose_name='Evidence', blank=True),
        ),
        migrations.AddField(
            model_name='entrychardata',
            name='evidence',
            field=cyano.history.HistoryManyToManyField(to='cyano.Evidence', verbose_name='Evidence', blank=True),
        ),
        migrations.AddField(
            model_name='entrybooleandata',
            name='evidence',
            field=cyano.history.HistoryManyToManyField(to='cyano.Evidence', verbose_name='Evidence', blank=True),
        ),
        migrations.AddField(
            model_name='entry',
            name='model_type',
            field=cyano.history.HistoryForeignKey(to='cyano.TableMeta'),
        ),
        migrations.AddField(
            model_name='entry',
            name='synonyms',
            field=cyano.history.HistoryManyToManyField(related_name='entry', verbose_name='Synonyms', to='cyano.Synonym', blank=True),
        ),
        migrations.AddField(
            model_name='dnafootprint',
            name='evidence',
            field=cyano.history.HistoryManyToManyField(to='cyano.Evidence', verbose_name='Evidence', blank=True),
        ),
        migrations.AddField(
            model_name='disulfidebond',
            name='evidence',
            field=cyano.history.HistoryManyToManyField(to='cyano.Evidence', verbose_name='Evidence', blank=True),
        ),
        migrations.AlterUniqueTogether(
            name='crossreference',
            unique_together=set([('xid', 'source')]),
        ),
        migrations.AddField(
            model_name='coenzymeparticipant',
            name='evidence',
            field=cyano.history.HistoryManyToManyField(to='cyano.Evidence', verbose_name='Evidence', blank=True),
        ),
        migrations.AddField(
            model_name='codon',
            name='evidence',
            field=cyano.history.HistoryManyToManyField(to='cyano.Evidence', verbose_name='Evidence', blank=True),
        ),
        migrations.AddField(
            model_name='biomasscomposition',
            name='evidence',
            field=cyano.history.HistoryManyToManyField(to='cyano.Evidence', verbose_name='Evidence', blank=True),
        ),
        migrations.AddField(
            model_name='bindingsite',
            name='evidence',
            field=cyano.history.HistoryManyToManyField(to='cyano.Evidence', verbose_name='Evidence', blank=True),
        ),
        migrations.AddField(
            model_name='basket',
            name='user',
            field=cyano.history.HistoryForeignKey(related_name='baskets', verbose_name='Users baskets', to='cyano.UserProfile'),
        ),
        migrations.CreateModel(
            name='ChromosomeFeature',
            fields=[
                ('parent_ptr_species_component', models.OneToOneField(parent_link=True, related_name='child_ptr_chromosome_feature', primary_key=True, serialize=False, to='cyano.SpeciesComponent', verbose_name='Species component')),
            ],
            options={
                'wid_unique': False,
                'concrete_entry_model': True,
                'facet_fields': ['type'],
                'verbose_name_plural': 'Chromosome features',
                'verbose_name': 'Chromosome feature',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'comments', 'publication_references'],
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), {'inline': 'positions'}, ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
            },
            bases=('cyano.speciescomponent',),
        ),
        migrations.CreateModel(
            name='Compartment',
            fields=[
                ('parent_ptr_species_component', models.OneToOneField(parent_link=True, related_name='child_ptr_compartment', primary_key=True, serialize=False, to='cyano.SpeciesComponent', verbose_name='Species component')),
            ],
            options={
                'wid_unique': False,
                'concrete_entry_model': True,
                'facet_fields': ['type'],
                'verbose_name_plural': 'Compartments',
                'verbose_name': 'Compartment',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'comments', 'publication_references'],
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Content', {'fields': [{'verbose_name': 'Metabolites (mM)', 'name': 'biomass_compositions'}, {'verbose_name': 'Protein monomers', 'name': 'protein_monomers'}, {'verbose_name': 'Protein complexes', 'name': 'protein_complexes'}]}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
            },
            bases=('cyano.speciescomponent',),
        ),
        migrations.CreateModel(
            name='MassSpectrometryJob',
            fields=[
                ('parent_ptr_species_component', models.OneToOneField(parent_link=True, related_name='child_ptr_mass_spectrometry', primary_key=True, serialize=False, to='cyano.SpeciesComponent', verbose_name='Species component')),
            ],
            options={
                'wid_unique': False,
                'concrete_entry_model': True,
                'facet_fields': ['type'],
                'verbose_name_plural': 'Mass Spectrometry Jobs',
                'verbose_name': 'Mass Spectrometry Job',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'comments', 'publication_references'],
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Mass Spectrometry', {'fields': [{'verbose_name': 'Target Peptides', 'name': 'target_peptide'}, {'verbose_name': 'Decoy Peptides', 'name': 'decoy_peptide'}, {'verbose_name': 'Related Proteins', 'name': 'related_proteins'}]}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
            },
            bases=('cyano.speciescomponent',),
        ),
        migrations.CreateModel(
            name='Molecule',
            fields=[
                ('parent_ptr_species_component', models.OneToOneField(parent_link=True, related_name='child_ptr_molecule', primary_key=True, serialize=False, to='cyano.SpeciesComponent', verbose_name='Species component')),
            ],
            options={
                'wid_unique': False,
                'concrete_entry_model': False,
                'facet_fields': ['type'],
                'verbose_name_plural': 'Molecules',
                'verbose_name': 'Molecule',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'comments', 'publication_references'],
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Function', {'fields': [{'verbose_name': 'Reaction participant', 'name': 'reaction_stoichiometry_participants'}, {'verbose_name': 'Complex subunit', 'name': 'protein_complex_biosythesis_participants'}]}), ('Parameters', {'fields': ['parameters']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
            },
            bases=('cyano.speciescomponent',),
        ),
        migrations.CreateModel(
            name='Note',
            fields=[
                ('parent_ptr_species_component', models.OneToOneField(parent_link=True, related_name='child_ptr_note', primary_key=True, serialize=False, to='cyano.SpeciesComponent', verbose_name='Species component')),
            ],
            options={
                'wid_unique': False,
                'concrete_entry_model': True,
                'facet_fields': ['type'],
                'verbose_name_plural': 'Notes',
                'verbose_name': 'Note',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'comments', 'publication_references'],
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
            },
            bases=('cyano.speciescomponent',),
        ),
        migrations.CreateModel(
            name='Parameter',
            fields=[
                ('parent_ptr_species_component', models.OneToOneField(parent_link=True, related_name='child_ptr_parameter', primary_key=True, serialize=False, to='cyano.SpeciesComponent', verbose_name='Species component')),
            ],
            options={
                'wid_unique': False,
                'concrete_entry_model': True,
                'facet_fields': ['type', 'reactions', 'molecules', 'state', 'process'],
                'verbose_name_plural': 'Misc. parameters',
                'verbose_name': 'Misc. parameter',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'value', 'reactions', 'molecules', 'state', 'process', 'comments', 'publication_references'],
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Value', {'fields': ['value']}), ('Associations', {'fields': ['reactions', 'molecules', 'state', 'process']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
            },
            bases=('cyano.speciescomponent',),
        ),
        migrations.CreateModel(
            name='Pathway',
            fields=[
                ('parent_ptr_species_component', models.OneToOneField(parent_link=True, related_name='child_ptr_pathway', primary_key=True, serialize=False, to='cyano.SpeciesComponent', verbose_name='Species component')),
            ],
            options={
                'wid_unique': True,
                'concrete_entry_model': True,
                'facet_fields': ['type'],
                'verbose_name_plural': 'Pathways',
                'group_field': 'type',
                'verbose_name': 'Pathway',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'comments', 'publication_references', 'type'],
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Navigator', {'fields': [{'verbose_name': 'Navigate to', 'name': 'navigator'}]}), ('Reactions', {'fields': [{'verbose_name': 'Reactions', 'name': 'reaction_map'}, 'reactions']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
            },
            bases=('cyano.speciescomponent',),
        ),
        migrations.CreateModel(
            name='Process',
            fields=[
                ('parent_ptr_species_component', models.OneToOneField(parent_link=True, related_name='child_ptr_process', primary_key=True, serialize=False, to='cyano.SpeciesComponent', verbose_name='Species component')),
                ('initialization_order', models.PositiveIntegerField(verbose_name='Initialization order')),
            ],
            options={
                'wid_unique': False,
                'concrete_entry_model': True,
                'facet_fields': ['type'],
                'verbose_name_plural': 'Processes',
                'verbose_name': 'Process',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'initialization_order', 'comments', 'publication_references'],
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Implementation', {'fields': ['initialization_order']}), ('Reactions', {'fields': [{'verbose_name': 'Chemical reactions', 'name': 'reactions'}, {'verbose_name': 'Complex formation reactions', 'name': 'formed_complexes'}]}), ('Parameters', {'fields': ['parameters']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
            },
            bases=('cyano.speciescomponent',),
        ),
        migrations.CreateModel(
            name='PublicationReference',
            fields=[
                ('parent_ptr_species_component', models.OneToOneField(parent_link=True, related_name='child_ptr_publicationreference', primary_key=True, serialize=False, to='cyano.SpeciesComponent', verbose_name='Species component')),
                ('authors', models.TextField(default='', verbose_name='Author(s)', blank=True)),
                ('editors', models.TextField(default='', verbose_name='Editor(s)', blank=True)),
                ('year', models.PositiveIntegerField(null=True, verbose_name='Year', blank=True)),
                ('title', models.TextField(default='', verbose_name='Title', blank=True)),
                ('publication', models.CharField(default='', max_length=255, verbose_name='Publication', blank=True)),
                ('publisher', models.CharField(default='', max_length=255, verbose_name='Publisher', blank=True)),
                ('volume', models.CharField(default='', max_length=255, verbose_name='Volume', blank=True)),
                ('issue', models.CharField(default='', max_length=255, verbose_name='Issue', blank=True)),
                ('pages', models.CharField(default='', max_length=255, verbose_name='Page(s)', blank=True)),
            ],
            options={
                'wid_unique': True,
                'concrete_entry_model': True,
                'facet_fields': ['type', 'year', 'publication'],
                'verbose_name_plural': 'Publication References',
                'verbose_name': 'Publication Reference',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'authors', 'editors', 'year', 'title', 'publication', 'publisher', 'volume', 'issue', 'pages', 'comments', 'publication_references'],
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Citation', {'fields': [{'verbose_name': 'Citation', 'name': 'citation'}]}), ('Cited by', {'fields': [{'verbose_name': 'Cited by', 'name': 'referenced_entries'}]}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
            },
            bases=('cyano.speciescomponent',),
        ),
        migrations.CreateModel(
            name='Reaction',
            fields=[
                ('parent_ptr_species_component', models.OneToOneField(parent_link=True, related_name='child_ptr_reaction', primary_key=True, serialize=False, to='cyano.SpeciesComponent', verbose_name='Species component')),
                ('direction', models.CharField(max_length=1, verbose_name='Direction', choices=[('f', 'Forward'), ('b', 'Backward'), ('r', 'Reversible')])),
                ('is_spontaneous', models.BooleanField(verbose_name='Is spontaneous (pH 7.5, 25C, <i>I</i> = 0)')),
                ('delta_g', models.FloatField(null=True, verbose_name='&Delta;G (pH 7.5, 25C, <i>I</i> = 0; kJ mol<sup>-1</sup>)', blank=True)),
                ('coenzymes', cyano.history.HistoryManyToManyField(related_name='reactions', verbose_name='Coenzymes', to='cyano.CoenzymeParticipant', blank=True)),
                ('enzyme', cyano.history.HistoryForeignKey(related_name='reactions', on_delete=django.db.models.deletion.SET_NULL, verbose_name='Enzyme', blank=True, to='cyano.EnzymeParticipant', null=True)),
                ('keq', cyano.history.HistoryForeignKey(related_name='+', on_delete=django.db.models.deletion.SET_NULL, verbose_name='K<sub>eq</sub>', blank=True, to='cyano.EntryPositiveFloatData', null=True)),
                ('kinetics_backward', cyano.history.HistoryForeignKey(related_name='reactions_backward', on_delete=django.db.models.deletion.SET_NULL, verbose_name='Backward kinetics', blank=True, to='cyano.Kinetics', null=True)),
                ('kinetics_forward', cyano.history.HistoryForeignKey(related_name='reactions_forward', on_delete=django.db.models.deletion.SET_NULL, verbose_name='Forward kinetics', blank=True, to='cyano.Kinetics', null=True)),
                ('map_coordinates', cyano.history.HistoryManyToManyField(related_name='reactions', verbose_name='Map coordinates', to='cyano.ReactionMapCoordinate', blank=True)),
                ('modification', cyano.history.HistoryForeignKey(related_name='reactions', on_delete=django.db.models.deletion.SET_NULL, verbose_name='Modification', blank=True, to='cyano.ModificationReaction', null=True)),
                ('optimal_ph', cyano.history.HistoryForeignKey(related_name='+', on_delete=django.db.models.deletion.SET_NULL, verbose_name='Optimal pH', blank=True, to='cyano.EntryPositiveFloatData', null=True)),
                ('optimal_temperature', cyano.history.HistoryForeignKey(related_name='+', on_delete=django.db.models.deletion.SET_NULL, verbose_name='Optimal temperature', blank=True, to='cyano.EntryFloatData', null=True)),
                ('pathways', cyano.history.HistoryManyToManyField(related_name='reactions', verbose_name='Pathways', to='cyano.Pathway', blank=True)),
                ('processes', cyano.history.HistoryForeignKey(related_name='reactions', on_delete=django.db.models.deletion.SET_NULL, verbose_name='Process', blank=True, to='cyano.Process', null=True)),
            ],
            options={
                'wid_unique': False,
                'concrete_entry_model': True,
                'facet_fields': ['type', 'direction', 'enzyme__protein', 'coenzymes__metabolite', 'is_spontaneous', 'pathways', 'processes', 'states'],
                'verbose_name_plural': 'Reactions',
                'verbose_name': 'Reaction',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'stoichiometry', 'direction', 'modification', 'enzyme', 'coenzymes', 'optimal_ph', 'optimal_temperature', 'is_spontaneous', 'delta_g', 'keq', 'kinetics_forward', 'kinetics_backward', 'pathways', 'processes', 'states', 'map_coordinates', 'comments', 'publication_references'],
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Reaction', {'fields': ['stoichiometry', 'modification']}), ('Catalysis', {'fields': ['enzyme', 'coenzymes', 'optimal_ph', 'optimal_temperature']}), ('Energetics', {'fields': ['is_spontaneous', 'delta_g', 'keq']}), ('Kinetics', {'fields': ['kinetics_forward', 'kinetics_backward']}), ('Parameters', {'fields': ['parameters']}), ('Associations', {'fields': ['pathways', 'processes', 'states']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
            },
            bases=('cyano.speciescomponent',),
        ),
        migrations.CreateModel(
            name='State',
            fields=[
                ('parent_ptr_species_component', models.OneToOneField(parent_link=True, related_name='child_ptr_state', primary_key=True, serialize=False, to='cyano.SpeciesComponent', verbose_name='Species component')),
            ],
            options={
                'wid_unique': False,
                'concrete_entry_model': True,
                'facet_fields': ['type'],
                'verbose_name_plural': 'States',
                'verbose_name': 'State',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'comments', 'publication_references'],
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Reactions', {'fields': ['reactions']}), ('Parameters', {'fields': ['parameters']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
            },
            bases=('cyano.speciescomponent',),
        ),
        migrations.CreateModel(
            name='TranscriptionalRegulation',
            fields=[
                ('parent_ptr_species_component', models.OneToOneField(parent_link=True, related_name='child_ptr_transcriptional_regulation', primary_key=True, serialize=False, to='cyano.SpeciesComponent', verbose_name='Species component')),
                ('activity', cyano.history.HistoryForeignKey(related_name='+', verbose_name='Fold-change activity', to='cyano.EntryPositiveFloatData')),
                ('affinity', cyano.history.HistoryForeignKey(related_name='+', on_delete=django.db.models.deletion.SET_NULL, verbose_name='Affinity', blank=True, to='cyano.EntryPositiveFloatData', null=True)),
                ('binding_site', cyano.history.HistoryForeignKey(related_name='transcriptional_regulations', on_delete=django.db.models.deletion.SET_NULL, verbose_name='Binding site', blank=True, to='cyano.BindingSite', null=True)),
            ],
            options={
                'wid_unique': False,
                'concrete_entry_model': True,
                'facet_fields': ['type', 'transcription_unit', 'transcription_factor'],
                'verbose_name_plural': 'Transcriptional regulation',
                'verbose_name': 'Transcriptional regulation',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'transcription_unit', 'transcription_factor', 'binding_site', 'affinity', 'activity', 'comments', 'publication_references'],
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Regulation', {'fields': ['transcription_unit', 'transcription_factor', 'binding_site', 'affinity', 'activity']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
            },
            bases=('cyano.speciescomponent',),
        ),
        migrations.CreateModel(
            name='Type',
            fields=[
                ('parent_ptr_species_component', models.OneToOneField(parent_link=True, related_name='child_ptr_type', primary_key=True, serialize=False, to='cyano.SpeciesComponent', verbose_name='Species component')),
            ],
            options={
                'wid_unique': True,
                'concrete_entry_model': True,
                'facet_fields': ['type', 'parent'],
                'verbose_name_plural': 'Types',
                'verbose_name': 'Type',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'parent', 'comments', 'publication_references'],
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type', 'parent', 'children', 'members']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
            },
            bases=('cyano.speciescomponent',),
        ),
        migrations.AddField(
            model_name='speciescomponent',
            name='cross_references',
            field=cyano.history.HistoryManyToManyField(related_name='cross_referenced_components', verbose_name='Cross references', to='cyano.CrossReference', blank=True),
        ),
        migrations.AddField(
            model_name='speciescomponent',
            name='parent',
            field=cyano.history.HistoryForeignKey(related_name='children', on_delete=django.db.models.deletion.SET_NULL, verbose_name='Parent', blank=True, to='cyano.SpeciesComponent', null=True),
        ),
        migrations.AddField(
            model_name='speciescomponent',
            name='species',
            field=cyano.history.HistoryForeignKey(related_name='components', verbose_name='Species', to='cyano.Species'),
        ),
        migrations.AddField(
            model_name='species',
            name='cross_references',
            field=cyano.history.HistoryManyToManyField(related_name='cross_referenced_species', verbose_name='Cross references', to='cyano.CrossReference', blank=True),
        ),
        migrations.AddField(
            model_name='evidence',
            name='species_component',
            field=cyano.history.HistoryManyToManyField(related_name='+', verbose_name='Species component', to='cyano.SpeciesComponent', blank=True),
        ),
        migrations.AlterUniqueTogether(
            name='entryuserobjectpermission',
            unique_together=set([('user', 'permission', 'content_object')]),
        ),
        migrations.AlterUniqueTogether(
            name='entrygroupobjectpermission',
            unique_together=set([('group', 'permission', 'content_object')]),
        ),
        migrations.AddField(
            model_name='basketcomponent',
            name='component',
            field=cyano.history.HistoryForeignKey(related_name='+', verbose_name='component', to='cyano.SpeciesComponent'),
        ),
        migrations.AddField(
            model_name='basketcomponent',
            name='species',
            field=cyano.history.HistoryForeignKey(related_name='+', verbose_name='Species component belongs to', to='cyano.Species'),
        ),
        migrations.CreateModel(
            name='Gene',
            fields=[
                ('parent_ptr_molecule', models.OneToOneField(parent_link=True, related_name='child_ptr_gene', primary_key=True, serialize=False, to='cyano.Molecule', verbose_name='Molecule')),
                ('symbol', models.CharField(default='', max_length=255, verbose_name='Symbol', blank=True)),
                ('coordinate', models.PositiveIntegerField(verbose_name='Coordinate (nt)')),
                ('length', models.PositiveIntegerField(verbose_name='Length (nt)')),
                ('direction', models.CharField(max_length=10, verbose_name='Direction', choices=[('f', 'Forward'), ('r', 'Reverse')])),
            ],
            options={
                'wid_unique': False,
                'concrete_entry_model': True,
                'facet_fields': ['type', 'chromosome', 'direction', 'is_essential', 'amino_acid'],
                'verbose_name_plural': 'Genes',
                'listing': ['wid', 'symbol'],
                'verbose_name': 'Gene',
                'field_list': ['id', 'wid', 'name', 'symbol', 'synonyms', 'cross_references', 'homologs', 'type', 'chromosome', 'coordinate', 'length', 'direction', 'is_essential', 'expression', 'half_life', 'codons', 'amino_acid', 'comments', 'publication_references'],
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'symbol', 'synonyms', {'verbose_name': 'Protein product', 'name': 'protein_monomers'}, 'cross_references', 'homologs']}), ('Classification', {'fields': ['type']}), ('Structure', {'fields': [{'verbose_name': 'Structure', 'name': 'structure'}, {'verbose_name': 'Structure Filter', 'name': 'structure_filter'}, {'verbose_name': 'Sequence', 'name': 'sequence'}, {'verbose_name': 'Transcription unit', 'name': 'transcription_units'}, {'verbose_name': 'Empirical formula (pH 7.5)', 'name': 'empirical_formula'}, {'verbose_name': 'Molecular weight (pH 7.5; Da)', 'name': 'molecular_weight'}]}), ('Functional genomics', {'fields': ['is_essential', 'expression', 'half_life', 'codons', 'amino_acid', {'verbose_name': 'Extinction coefficient <br/>(260 nm, 25C, pH 7.0)', 'name': 'extinction_coefficient'}, {'verbose_name': 'pI', 'name': 'pi'}]}), ('Function', {'fields': [{'verbose_name': 'Reaction participant', 'name': 'reaction_stoichiometry_participants'}, {'verbose_name': 'Complex subunit', 'name': 'protein_complex_biosythesis_participants'}]}), ('Parameters', {'fields': ['parameters']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
            },
            bases=('cyano.molecule',),
        ),
        migrations.CreateModel(
            name='Genome',
            fields=[
                ('parent_ptr_molecule', models.OneToOneField(parent_link=True, related_name='child_ptr_genome', primary_key=True, serialize=False, to='cyano.Molecule', verbose_name='Molecule')),
                ('sequence', models.TextField(default='', blank=True, verbose_name='Sequence', validators=[cyano.models.validate_dna_sequence])),
                ('length', models.PositiveIntegerField(verbose_name='Length (nt)')),
            ],
            options={
                'wid_unique': False,
                'concrete_entry_model': True,
                'facet_fields': ['type'],
                'verbose_name_plural': 'Genome',
                'verbose_name': 'Genome',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'sequence', 'length', 'comments', 'publication_references'],
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Sequence', {'fields': [{'verbose_name': 'Structure Filter', 'name': 'structure_filter'}, {'verbose_name': 'Structure', 'name': 'structure'}, {'verbose_name': 'Extinction coefficient <br/>(260 nm, 25C, pH 7.0)', 'name': 'extinction_coefficient'}, {'verbose_name': 'pI', 'name': 'pi'}]}), ('Features', {'fields': ['genes', {'verbose_name': 'Other features', 'name': 'features'}]}), ('Function', {'fields': [{'verbose_name': 'Reaction participant', 'name': 'reaction_stoichiometry_participants'}, {'verbose_name': 'Complex subunit', 'name': 'protein_complex_biosythesis_participants'}]}), ('Parameters', {'fields': ['parameters']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
            },
            bases=('cyano.molecule',),
        ),
        migrations.CreateModel(
            name='Metabolite',
            fields=[
                ('parent_ptr_molecule', models.OneToOneField(parent_link=True, related_name='child_ptr_metabolite', primary_key=True, serialize=False, to='cyano.Molecule', verbose_name='Molecule')),
                ('traditional_name', models.TextField(default='', verbose_name='Traditional name', blank=True)),
                ('iupac_name', models.TextField(default='', verbose_name='IUPAC name', blank=True)),
                ('empirical_formula', models.TextField(verbose_name='Empirical formula (pH 7.5)', validators=[django.core.validators.RegexValidator(regex='^([A-Z][a-z]*[0-9]*)+$', message='Invalid empirical formula')])),
                ('smiles', models.TextField(default='', verbose_name='SMILES (pH 7.5)', blank=True)),
                ('charge', models.IntegerField(verbose_name='Charge (pH 7.5)')),
                ('is_hydrophobic', models.BooleanField(verbose_name='Is hydrophobic')),
                ('volume', models.FloatField(blank=True, null=True, verbose_name='van der Waals volume <br/>(pH 7.5; &#8491;<sup>3</sup> molecule<sup>-1</sup>)', validators=[django.core.validators.MinValueValidator(0)])),
                ('deltag_formation', models.FloatField(null=True, verbose_name="&Delta;<sub>f</sub>G<sup>'o</sup> (pH 7.5, 25C, I = 0; kJ mol<sup>-1</sup>)", blank=True)),
                ('pka', models.FloatField(blank=True, null=True, verbose_name='pK<sub>a</sub>', validators=[django.core.validators.MinValueValidator(0)])),
                ('pi', models.FloatField(blank=True, null=True, verbose_name='pI', validators=[django.core.validators.MinValueValidator(0)])),
                ('log_p', models.FloatField(null=True, verbose_name='logP', blank=True)),
                ('log_d', models.FloatField(null=True, verbose_name='logD (pH 7.5)', blank=True)),
            ],
            options={
                'wid_unique': False,
                'concrete_entry_model': True,
                'facet_fields': ['type', 'charge', 'is_hydrophobic'],
                'verbose_name_plural': 'Metabolites',
                'verbose_name': 'Metabolite',
                'field_list': ['id', 'wid', 'name', 'traditional_name', 'iupac_name', 'synonyms', 'cross_references', 'type', 'empirical_formula', 'smiles', 'charge', 'is_hydrophobic', 'volume', 'deltag_formation', 'pka', 'pi', 'log_p', 'log_d', 'biomass_composition', 'media_composition', 'map_coordinates', 'comments', 'publication_references'],
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'traditional_name', 'iupac_name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Structure', {'fields': [{'verbose_name': 'Structure', 'name': 'structure'}, 'empirical_formula', 'smiles', 'charge', 'is_hydrophobic', {'verbose_name': 'Molecular weight (Da)', 'name': 'molecular_weight'}, 'volume', 'deltag_formation', 'pka', 'pi', 'log_p', 'log_d']}), ('Concentrations', {'fields': ['biomass_composition', 'media_composition']}), ('Function', {'fields': [{'verbose_name': 'Coenzyme', 'name': 'coenzyme_participants'}, {'verbose_name': 'Complex subunit', 'name': 'protein_complex_biosythesis_participants'}, {'verbose_name': 'Prosthetic group', 'name': 'prosthetic_group_participants'}, {'verbose_name': 'Reaction participant', 'name': 'reaction_stoichiometry_participants'}]}), ('Parameters', {'fields': ['parameters']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
            },
            bases=('cyano.molecule',),
        ),
        migrations.CreateModel(
            name='Protein',
            fields=[
                ('parent_ptr_molecule', models.OneToOneField(parent_link=True, related_name='child_ptr_protein', primary_key=True, serialize=False, to='cyano.Molecule', verbose_name='Molecule')),
            ],
            options={
                'wid_unique': False,
                'concrete_entry_model': False,
                'facet_fields': ['type', 'chaperones', 'dna_footprint__binding', 'dna_footprint__region'],
                'verbose_name_plural': 'Proteins',
                'verbose_name': 'Protein',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'prosthetic_groups', 'chaperones', 'dna_footprint', 'regulatory_rule', 'comments', 'publication_references'],
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Structure', {'fields': ['prosthetic_groups', 'chaperones', 'dna_footprint']}), ('Regulation', {'fields': ['regulatory_rule']}), ('Function', {'fields': [{'verbose_name': 'Enzyme', 'name': 'enzyme_participants'}, {'verbose_name': 'Transcriptional regulation', 'name': 'transcriptional_regulations'}, {'verbose_name': 'Protein folding substrates', 'name': 'chaperone_substrates'}, {'verbose_name': 'Reaction participant', 'name': 'reaction_stoichiometry_participants'}, {'verbose_name': 'Complex subunit', 'name': 'protein_complex_biosythesis_participants'}]}), ('Parameters', {'fields': ['parameters']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
            },
            bases=('cyano.molecule',),
        ),
        migrations.CreateModel(
            name='Stimulus',
            fields=[
                ('parent_ptr_molecule', models.OneToOneField(parent_link=True, related_name='child_ptr_stimulus', primary_key=True, serialize=False, to='cyano.Molecule', verbose_name='Molecule')),
                ('value', cyano.history.HistoryForeignKey(related_name='+', verbose_name='Value', to='cyano.EntryFloatData')),
            ],
            options={
                'wid_unique': False,
                'concrete_entry_model': True,
                'facet_fields': ['type', 'value__units'],
                'verbose_name_plural': 'Stimuli',
                'verbose_name': 'Stimulus',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'value', 'comments', 'publication_references'],
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Value', {'fields': ['value']}), ('Function', {'fields': [{'verbose_name': 'Reaction participant', 'name': 'reaction_stoichiometry_participants'}, {'verbose_name': 'Complex subunit', 'name': 'protein_complex_biosythesis_participants'}]}), ('Parameters', {'fields': ['parameters']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
            },
            bases=('cyano.molecule',),
        ),
        migrations.CreateModel(
            name='TranscriptionUnit',
            fields=[
                ('parent_ptr_molecule', models.OneToOneField(parent_link=True, related_name='child_ptr_transcription_unit', primary_key=True, serialize=False, to='cyano.Molecule', verbose_name='Molecule')),
                ('promoter_35_coordinate', models.IntegerField(null=True, verbose_name='Promoter -35 box coordinate (nt)', blank=True)),
                ('promoter_35_length', models.IntegerField(null=True, verbose_name='Promoter -35 box length (nt)', blank=True)),
                ('promoter_10_coordinate', models.IntegerField(null=True, verbose_name='Promoter -10 box coordinate (nt)', blank=True)),
                ('promoter_10_length', models.IntegerField(null=True, verbose_name='Promoter -10 box length (nt)', blank=True)),
                ('tss_coordinate', models.IntegerField(null=True, verbose_name='Transcription start site coordinate (nt)', blank=True)),
                ('genes', cyano.history.HistoryManyToManyField(related_name='transcription_units', verbose_name='Genes', to='cyano.Gene')),
            ],
            options={
                'wid_unique': False,
                'concrete_entry_model': True,
                'facet_fields': ['type'],
                'verbose_name_plural': 'Transcription units',
                'verbose_name': 'Transcription unit',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'genes', 'promoter_35_coordinate', 'promoter_35_length', 'promoter_10_coordinate', 'promoter_10_length', 'tss_coordinate', 'comments', 'publication_references'],
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Structure (Hayflick media, 37C)', {'fields': [{'verbose_name': 'Structure', 'name': 'structure'}, {'verbose_name': 'Structure Filter', 'name': 'structure_filter'}, 'genes', 'promoter_35_coordinate', 'promoter_35_length', 'promoter_10_coordinate', 'promoter_10_length', 'tss_coordinate', {'verbose_name': 'Sequence', 'name': 'sequence'}]}), ('Regulation', {'fields': [{'verbose_name': 'Regulation', 'name': 'transcriptional_regulations'}]}), ('Function', {'fields': [{'verbose_name': 'Reaction participant', 'name': 'reaction_stoichiometry_participants'}, {'verbose_name': 'Complex subunit', 'name': 'protein_complex_biosythesis_participants'}]}), ('Parameters', {'fields': ['parameters']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
            },
            bases=('cyano.molecule',),
        ),
        migrations.AddField(
            model_name='speciescomponent',
            name='publication_references',
            field=cyano.history.HistoryManyToManyField(related_name='publication_referenced_components', verbose_name='Publications', to='cyano.PublicationReference', blank=True),
        ),
        migrations.AddField(
            model_name='speciescomponent',
            name='type',
            field=cyano.history.HistoryManyToManyField(related_name='members', verbose_name='Type', to='cyano.Type', blank=True),
        ),
        migrations.AddField(
            model_name='species',
            name='publication_references',
            field=cyano.history.HistoryManyToManyField(related_name='publication_referenced_species', verbose_name='Publication references', to='cyano.PublicationReference', blank=True),
        ),
        migrations.AddField(
            model_name='reactionstoichiometryparticipant',
            name='compartment',
            field=cyano.history.HistoryForeignKey(related_name='+', verbose_name='Compartment', to='cyano.Compartment'),
        ),
        migrations.AddField(
            model_name='reactionstoichiometryparticipant',
            name='molecule',
            field=cyano.history.HistoryForeignKey(related_name='reaction_stoichiometry_participants', verbose_name='Molecule', to='cyano.Molecule'),
        ),
        migrations.AddField(
            model_name='reaction',
            name='states',
            field=cyano.history.HistoryForeignKey(related_name='reactions', on_delete=django.db.models.deletion.SET_NULL, verbose_name='State', blank=True, to='cyano.State', null=True),
        ),
        migrations.AddField(
            model_name='reaction',
            name='stoichiometry',
            field=cyano.history.HistoryManyToManyField(related_name='reactions', verbose_name='Stoichiometry', to='cyano.ReactionStoichiometryParticipant'),
        ),
        migrations.AddField(
            model_name='proteincomplexbiosythesisparticipant',
            name='compartment',
            field=cyano.history.HistoryForeignKey(related_name='+', verbose_name='Compartment', to='cyano.Compartment'),
        ),
        migrations.AddField(
            model_name='proteincomplexbiosythesisparticipant',
            name='molecule',
            field=cyano.history.HistoryForeignKey(related_name='protein_complex_biosythesis_participants', verbose_name='Molecule', to='cyano.Molecule'),
        ),
        migrations.AddField(
            model_name='prostheticgroupparticipant',
            name='compartment',
            field=cyano.history.HistoryForeignKey(related_name='+', verbose_name='Compartment', to='cyano.Compartment'),
        ),
        migrations.AddField(
            model_name='parameter',
            name='molecules',
            field=cyano.history.HistoryManyToManyField(related_name='parameters', verbose_name='Molecules', to='cyano.Molecule', blank=True),
        ),
        migrations.AddField(
            model_name='parameter',
            name='process',
            field=cyano.history.HistoryForeignKey(related_name='parameters', on_delete=django.db.models.deletion.SET_NULL, verbose_name='Process', blank=True, to='cyano.Process', null=True),
        ),
        migrations.AddField(
            model_name='parameter',
            name='reactions',
            field=cyano.history.HistoryManyToManyField(related_name='parameters', verbose_name='Reactions', to='cyano.Reaction', blank=True),
        ),
        migrations.AddField(
            model_name='parameter',
            name='state',
            field=cyano.history.HistoryForeignKey(related_name='parameters', on_delete=django.db.models.deletion.SET_NULL, verbose_name='State', blank=True, to='cyano.State', null=True),
        ),
        migrations.AddField(
            model_name='parameter',
            name='value',
            field=cyano.history.HistoryForeignKey(verbose_name='Value', to='cyano.EntryCharData'),
        ),
        migrations.AddField(
            model_name='modificationreaction',
            name='compartment',
            field=cyano.history.HistoryForeignKey(related_name='+', verbose_name='Compartment', to='cyano.Compartment'),
        ),
        migrations.AddField(
            model_name='modificationreaction',
            name='molecule',
            field=cyano.history.HistoryForeignKey(related_name='modification_reactions', verbose_name='Molecule', to='cyano.Molecule'),
        ),
        migrations.AddField(
            model_name='metabolitemapcoordinate',
            name='compartment',
            field=cyano.history.HistoryForeignKey(related_name='+', verbose_name='Compartment', to='cyano.Compartment'),
        ),
        migrations.AddField(
            model_name='featureposition',
            name='chromosome_feature',
            field=cyano.history.HistoryForeignKey(related_name='positions', verbose_name='', to='cyano.ChromosomeFeature'),
        ),
        migrations.AddField(
            model_name='evidence',
            name='references',
            field=cyano.history.HistoryForeignKey(related_name='evidence', verbose_name='References', blank=True, to='cyano.PublicationReference', null=True),
        ),
        migrations.AddField(
            model_name='enzymeparticipant',
            name='compartment',
            field=cyano.history.HistoryForeignKey(related_name='+', verbose_name='Compartment', to='cyano.Compartment'),
        ),
        migrations.AddField(
            model_name='coenzymeparticipant',
            name='compartment',
            field=cyano.history.HistoryForeignKey(related_name='+', verbose_name='Compartment', to='cyano.Compartment'),
        ),
        migrations.AddField(
            model_name='biomasscomposition',
            name='compartment',
            field=cyano.history.HistoryForeignKey(related_name='biomass_compositions', verbose_name='Compartment', to='cyano.Compartment'),
        ),
        migrations.CreateModel(
            name='Chromosome',
            fields=[
                ('parent_ptr_genome', models.OneToOneField(parent_link=True, related_name='child_ptr_chromosome', primary_key=True, serialize=False, to='cyano.Genome', verbose_name='Genome')),
            ],
            options={
                'wid_unique': False,
                'concrete_entry_model': True,
                'facet_fields': ['type'],
                'verbose_name_plural': 'Chromosomes',
                'verbose_name': 'Chromosome',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'sequence', 'length', 'comments', 'publication_references'],
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Sequence', {'fields': [{'verbose_name': 'Structure Filter', 'name': 'structure_filter'}, {'verbose_name': 'Structure', 'name': 'structure'}, {'verbose_name': 'Extinction coefficient <br/>(260 nm, 25C, pH 7.0)', 'name': 'extinction_coefficient'}, {'verbose_name': 'pI', 'name': 'pi'}]}), ('Features', {'fields': ['genes', {'verbose_name': 'Other features', 'name': 'features'}]}), ('Function', {'fields': [{'verbose_name': 'Reaction participant', 'name': 'reaction_stoichiometry_participants'}, {'verbose_name': 'Complex subunit', 'name': 'protein_complex_biosythesis_participants'}]}), ('Parameters', {'fields': ['parameters']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
            },
            bases=('cyano.genome',),
        ),
        migrations.CreateModel(
            name='MassSpectrometryProtein',
            fields=[
                ('protein_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='cyano.Protein')),
                ('parent_ptr_protein', models.OneToOneField(parent_link=True, related_name='child_ptr_ms_protein', verbose_name='Species component', to='cyano.SpeciesComponent')),
                ('score', models.FloatField(verbose_name='Protein Score')),
                ('coverage', models.FloatField(verbose_name='% Coverage')),
                ('sequence', models.TextField(verbose_name='Sequence', validators=[cyano.models.validate_protein_sequence])),
                ('length', models.PositiveIntegerField(verbose_name='Length (nt)')),
                ('pi', models.FloatField(verbose_name='Protein PI')),
                ('mass', models.FloatField(verbose_name='Protein Mass (Da)')),
                ('ambiguous', cyano.history.HistoryManyToManyField(related_name='ambiguous', verbose_name='Ambiguous Proteins', to='cyano.EntryBasicTextData')),
                ('sub', cyano.history.HistoryManyToManyField(related_name='sub', verbose_name='Sub-Proteins', to='cyano.EntryBasicTextData')),
            ],
            options={
                'wid_unique': False,
                'concrete_entry_model': True,
                'facet_fields': ['type', 'chaperones', 'dna_footprint__binding', 'dna_footprint__region'],
                'verbose_name_plural': 'Mass Spectrometry Proteins',
                'verbose_name': 'Mass Spectrometry Protein',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'prosthetic_groups', 'chaperones', 'dna_footprint', 'regulatory_rule', 'score', 'coverage', 'pi', 'mass', 'sequence', 'comments', 'publication_references'],
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Related', {'fields': ['parent', {'verbose_name': 'Ambiguous Proteins', 'name': 'ambiguous_proteins'}, {'verbose_name': 'Sun-Proteins', 'name': 'sub_proteins'}]}), ('Structure', {'fields': ['prosthetic_groups', 'chaperones', 'dna_footprint', {'verbose_name': 'Sequence', 'name': 'sequence'}, {'verbose_name': 'Structure', 'name': 'structure'}]}), ('Regulation', {'fields': ['regulatory_rule']}), ('Function', {'fields': [{'verbose_name': 'Enzyme', 'name': 'enzyme_participants'}, {'verbose_name': 'Transcriptional regulation', 'name': 'transcriptional_regulations'}, {'verbose_name': 'Protein folding substrates', 'name': 'chaperone_substrates'}, {'verbose_name': 'Reaction participant', 'name': 'reaction_stoichiometry_participants'}, {'verbose_name': 'Complex subunit', 'name': 'protein_complex_biosythesis_participants'}]}), ('Statistics', {'fields': ['score', 'coverage', 'pi', 'mass']}), {'inline': 'protein_details'}, ('Parameters', {'fields': ['parameters']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
            },
            bases=('cyano.protein',),
        ),
        migrations.CreateModel(
            name='Peptide',
            fields=[
                ('protein_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='cyano.Protein')),
                ('parent_ptr_protein', models.OneToOneField(parent_link=True, related_name='child_ptr_peptide', verbose_name='Species component', to='cyano.SpeciesComponent')),
                ('sequence', models.TextField(default='', blank=True, verbose_name='Sequence', validators=[cyano.models.validate_protein_sequence])),
                ('length', models.PositiveIntegerField(verbose_name='Length (nt)')),
                ('proteotypic', models.NullBooleanField(verbose_name='Proteotypic')),
                ('charge', models.IntegerField(verbose_name='Charge')),
                ('mass', models.FloatField(verbose_name='m/z')),
                ('zscore', models.FloatField(verbose_name='z-score')),
                ('retention_time', models.FloatField(verbose_name='Retention Time')),
                ('proteins', cyano.history.HistoryManyToManyField(related_name='peptides', verbose_name='Proteins belonging to the Peptide', to='cyano.EntryBasicTextData')),
            ],
            options={
                'wid_unique': False,
                'concrete_entry_model': True,
                'facet_fields': ['type', 'chaperones', 'dna_footprint__binding', 'dna_footprint__region'],
                'verbose_name_plural': 'Peptides',
                'verbose_name': 'Peptide',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'prosthetic_groups', 'chaperones', 'dna_footprint', 'regulatory_rule', 'comments', 'publication_references'],
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Related', {'fields': ['parent', {'verbose_name': 'Matched Proteins', 'name': 'matched_proteins'}]}), ('Structure', {'fields': ['prosthetic_groups', 'chaperones', 'dna_footprint', {'verbose_name': 'Sequence', 'name': 'sequence'}]}), ('Regulation', {'fields': ['regulatory_rule']}), ('Function', {'fields': [{'verbose_name': 'Enzyme', 'name': 'enzyme_participants'}, {'verbose_name': 'Transcriptional regulation', 'name': 'transcriptional_regulations'}, {'verbose_name': 'Protein folding substrates', 'name': 'chaperone_substrates'}, {'verbose_name': 'Reaction participant', 'name': 'reaction_stoichiometry_participants'}, {'verbose_name': 'Complex subunit', 'name': 'protein_complex_biosythesis_participants'}]}), ('Statistics', {'fields': ['proteotypic', 'charge', 'mass', 'zscore', 'retention_time']}), ('Parameters', {'fields': ['parameters']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
            },
            bases=('cyano.protein',),
        ),
        migrations.CreateModel(
            name='Plasmid',
            fields=[
                ('parent_ptr_genome', models.OneToOneField(parent_link=True, related_name='child_ptr_plasmid', primary_key=True, serialize=False, to='cyano.Genome', verbose_name='Genome')),
            ],
            options={
                'wid_unique': False,
                'concrete_entry_model': True,
                'facet_fields': ['type'],
                'verbose_name_plural': 'Plasmids',
                'verbose_name': 'Plasmid',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'sequence', 'length', 'comments', 'publication_references'],
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Sequence', {'fields': [{'verbose_name': 'Structure Filter', 'name': 'structure_filter'}, {'verbose_name': 'Structure', 'name': 'structure'}, {'verbose_name': 'Extinction coefficient <br/>(260 nm, 25C, pH 7.0)', 'name': 'extinction_coefficient'}, {'verbose_name': 'pI', 'name': 'pi'}]}), ('Features', {'fields': ['genes', {'verbose_name': 'Other features', 'name': 'features'}]}), ('Function', {'fields': [{'verbose_name': 'Reaction participant', 'name': 'reaction_stoichiometry_participants'}, {'verbose_name': 'Complex subunit', 'name': 'protein_complex_biosythesis_participants'}]}), ('Parameters', {'fields': ['parameters']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
            },
            bases=('cyano.genome',),
        ),
        migrations.CreateModel(
            name='ProteinComplex',
            fields=[
                ('parent_ptr_protein', models.OneToOneField(parent_link=True, related_name='child_ptr_protein_complex', primary_key=True, serialize=False, to='cyano.Protein', verbose_name='Protein')),
                ('biosynthesis', cyano.history.HistoryManyToManyField(related_name='protein_complexes', verbose_name='Biosynthesis', to='cyano.ProteinComplexBiosythesisParticipant')),
                ('disulfide_bonds', cyano.history.HistoryManyToManyField(related_name='protein_complexes', verbose_name='Disulfide bonds (pH 7.5)', to='cyano.DisulfideBond', blank=True)),
                ('formation_process', cyano.history.HistoryForeignKey(related_name='formed_complexes', on_delete=django.db.models.deletion.SET_NULL, verbose_name='Formation process', blank=True, to='cyano.Process', null=True)),
            ],
            options={
                'wid_unique': False,
                'concrete_entry_model': True,
                'facet_fields': ['type', 'dna_footprint__binding', 'dna_footprint__region', 'formation_process', 'chaperones'],
                'verbose_name_plural': 'Protein complexes',
                'verbose_name': 'Protein complex',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'biosynthesis', 'disulfide_bonds', 'prosthetic_groups', 'dna_footprint', 'formation_process', 'chaperones', 'regulatory_rule', 'comments', 'publication_references'],
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Structure', {'fields': ['biosynthesis', {'verbose_name': 'No. subunits', 'name': 'num_subunits'}, 'disulfide_bonds', 'prosthetic_groups', 'dna_footprint', {'verbose_name': 'Empirical formula (pH 7.5)', 'name': 'empirical_formula'}, {'verbose_name': 'Molecular weight (pH 7.5; Da)', 'name': 'molecular_weight'}, {'verbose_name': 'Extinction coefficient <br/>(260 nm, 25C, pH 7.0)', 'name': 'extinction_coefficient'}, {'verbose_name': 'Half life (OD (600 nm) = 0.3, <br/>M9 media, 36C; min)', 'name': 'half_life'}]}), ('Synthesis', {'fields': ['formation_process', 'chaperones', {'verbose_name': 'Localization', 'name': 'localization'}]}), ('Regulation', {'fields': ['regulatory_rule']}), ('Function', {'fields': [{'verbose_name': 'Enzyme', 'name': 'enzyme_participants'}, {'verbose_name': 'Transcriptional regulation', 'name': 'transcriptional_regulations'}, {'verbose_name': 'Protein folding substrates', 'name': 'chaperone_substrates'}, {'verbose_name': 'Reaction participant', 'name': 'reaction_stoichiometry_participants'}, {'verbose_name': 'Complex subunit', 'name': 'protein_complex_biosythesis_participants'}]}), ('Parameters', {'fields': ['parameters']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
            },
            bases=('cyano.protein',),
        ),
        migrations.CreateModel(
            name='ProteinMonomer',
            fields=[
                ('parent_ptr_protein', models.OneToOneField(parent_link=True, related_name='child_ptr_protein_monomer', primary_key=True, serialize=False, to='cyano.Protein', verbose_name='Protein')),
            ],
            options={
                'wid_unique': False,
                'concrete_entry_model': True,
                'facet_fields': ['type', 'is_n_terminal_methionine_cleaved__value', 'signal_sequence__type', 'signal_sequence__location', 'dna_footprint__binding', 'dna_footprint__region', 'localization', 'chaperones'],
                'verbose_name_plural': 'Protein monomers',
                'verbose_name': 'Protein monomer',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'gene', 'is_n_terminal_methionine_cleaved', 'signal_sequence', 'prosthetic_groups', 'dna_footprint', 'localization', 'chaperones', 'regulatory_rule', 'comments', 'publication_references'],
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Genetics', {'fields': ['gene']}), ('Structure', {'fields': [{'verbose_name': 'Sequence', 'name': 'sequence'}, 'is_n_terminal_methionine_cleaved', 'signal_sequence', 'prosthetic_groups', 'disulfide_bonds', 'dna_footprint', {'verbose_name': 'Empirical formula (pH 7.5)', 'name': 'empirical_formula'}, {'verbose_name': 'Molecular weight (pH 7.5; Da)', 'name': 'molecular_weight'}, {'verbose_name': 'Extinction coefficient <br/>(260 nm, 25C, pH 7.0)', 'name': 'extinction_coefficient'}, {'verbose_name': 'pI', 'name': 'pi'}, {'verbose_name': 'Instability index', 'name': 'instability'}, {'verbose_name': 'Is stable', 'name': 'is_stable'}, {'verbose_name': 'Aliphatic index', 'name': 'aliphatic'}, {'verbose_name': 'GRAVY (25C, pH 7.0)', 'name': 'gravy'}, {'verbose_name': 'Half life (OD (600 nm) = 0.3, <br/>M9 media, 36C; min)', 'name': 'half_life'}]}), ('Synthesis', {'fields': ['localization', 'chaperones']}), ('Regulation', {'fields': ['regulatory_rule']}), ('Function', {'fields': [{'verbose_name': 'Enzyme', 'name': 'enzyme_participants'}, {'verbose_name': 'Transcriptional regulation', 'name': 'transcriptional_regulations'}, {'verbose_name': 'Protein folding substrates', 'name': 'chaperone_substrates'}, {'verbose_name': 'Reaction participant', 'name': 'reaction_stoichiometry_participants'}, {'verbose_name': 'Complex subunit', 'name': 'protein_complex_biosythesis_participants'}]}), ('Parameters', {'fields': ['parameters']}), ('Interactions', {'fields': [{'verbose_name': 'Protein/Metabolite interactions', 'name': 'interactions'}]}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
            },
            bases=('cyano.protein',),
        ),
        migrations.AddField(
            model_name='transcriptionalregulation',
            name='transcription_factor',
            field=cyano.history.HistoryForeignKey(related_name='transcriptional_regulations', verbose_name='Transcripton factor', to='cyano.Protein'),
        ),
        migrations.AddField(
            model_name='transcriptionalregulation',
            name='transcription_unit',
            field=cyano.history.HistoryForeignKey(related_name='transcriptional_regulations', verbose_name='Transcription unit', to='cyano.TranscriptionUnit'),
        ),
        migrations.AddField(
            model_name='protein',
            name='chaperones',
            field=cyano.history.HistoryManyToManyField(related_name='chaperone_substrates', verbose_name='Chaperones', to='cyano.Protein', blank=True),
        ),
        migrations.AddField(
            model_name='protein',
            name='dna_footprint',
            field=cyano.history.HistoryForeignKey(related_name='proteins', verbose_name='DNA footprint', blank=True, to='cyano.DNAFootprint', null=True),
        ),
        migrations.AddField(
            model_name='protein',
            name='prosthetic_groups',
            field=cyano.history.HistoryManyToManyField(related_name='proteins', verbose_name='Prosthetic groups', to='cyano.ProstheticGroupParticipant', blank=True),
        ),
        migrations.AddField(
            model_name='protein',
            name='regulatory_rule',
            field=cyano.history.HistoryForeignKey(related_name='+', on_delete=django.db.models.deletion.SET_NULL, verbose_name='Regulatory rule', blank=True, to='cyano.EntryCharData', null=True),
        ),
        migrations.AddField(
            model_name='prostheticgroupparticipant',
            name='metabolite',
            field=cyano.history.HistoryForeignKey(related_name='prosthetic_group_participants', verbose_name='Metabolite', to='cyano.Metabolite'),
        ),
        migrations.AddField(
            model_name='metabolite',
            name='biomass_composition',
            field=cyano.history.HistoryManyToManyField(related_name='metabolites', verbose_name='Biomass composition (SP4 media, <br/>5% CO<sub>2</sub>, 37C; mmol gDCW<sup>-1</sup>)', to='cyano.BiomassComposition', blank=True),
        ),
        migrations.AddField(
            model_name='metabolite',
            name='map_coordinates',
            field=cyano.history.HistoryManyToManyField(related_name='metabolites', verbose_name='Map coordinates', to='cyano.MetaboliteMapCoordinate', blank=True),
        ),
        migrations.AddField(
            model_name='metabolite',
            name='media_composition',
            field=cyano.history.HistoryForeignKey(related_name='metabolites', on_delete=django.db.models.deletion.SET_NULL, verbose_name='Media composition (SP4; mM)', blank=True, to='cyano.MediaComposition', null=True),
        ),
        migrations.AddField(
            model_name='gene',
            name='amino_acid',
            field=cyano.history.HistoryForeignKey(related_name='genes', on_delete=django.db.models.deletion.SET_NULL, verbose_name='Amino acid', blank=True, to='cyano.Metabolite', null=True),
        ),
        migrations.AddField(
            model_name='gene',
            name='chromosome',
            field=cyano.history.HistoryForeignKey(related_name='genes', verbose_name='Chromosome or Plasmid', to='cyano.Genome'),
        ),
        migrations.AddField(
            model_name='gene',
            name='codons',
            field=cyano.history.HistoryManyToManyField(related_name='genes', verbose_name='Codons', to='cyano.Codon', blank=True),
        ),
        migrations.AddField(
            model_name='gene',
            name='expression',
            field=cyano.history.HistoryForeignKey(related_name='+', verbose_name='Relative expression', blank=True, to='cyano.EntryPositiveFloatData', null=True),
        ),
        migrations.AddField(
            model_name='gene',
            name='half_life',
            field=cyano.history.HistoryForeignKey(related_name='+', verbose_name='Half life', blank=True, to='cyano.EntryPositiveFloatData', null=True),
        ),
        migrations.AddField(
            model_name='gene',
            name='homologs',
            field=cyano.history.HistoryManyToManyField(related_name='genes', verbose_name='Homologs', to='cyano.Homolog', blank=True),
        ),
        migrations.AddField(
            model_name='gene',
            name='is_essential',
            field=cyano.history.HistoryForeignKey(related_name='+', verbose_name='Is essential', blank=True, to='cyano.EntryBooleanData', null=True),
        ),
        migrations.AddField(
            model_name='featureposition',
            name='chromosome',
            field=cyano.history.HistoryForeignKey(related_name='features', verbose_name='Chromosome or Plasmid', to='cyano.Genome'),
        ),
        migrations.AddField(
            model_name='enzymeparticipant',
            name='protein',
            field=cyano.history.HistoryForeignKey(related_name='enzyme_participants', verbose_name='Protein', to='cyano.Protein'),
        ),
        migrations.AddField(
            model_name='coenzymeparticipant',
            name='metabolite',
            field=cyano.history.HistoryForeignKey(related_name='coenzyme_participants', verbose_name='Metabolite', to='cyano.Metabolite'),
        ),
        migrations.AddField(
            model_name='proteinmonomer',
            name='gene',
            field=cyano.history.HistoryForeignKey(related_name='protein_monomers', verbose_name='Gene', to='cyano.Gene'),
        ),
        migrations.AddField(
            model_name='proteinmonomer',
            name='is_n_terminal_methionine_cleaved',
            field=cyano.history.HistoryForeignKey(related_name='+', verbose_name='Is N-terminal methionine cleaved', to='cyano.EntryBooleanData', null=True),
        ),
        migrations.AddField(
            model_name='proteinmonomer',
            name='localization',
            field=cyano.history.HistoryForeignKey(related_name='protein_monomers', verbose_name='Localization', to='cyano.Compartment', null=True),
        ),
        migrations.AddField(
            model_name='proteinmonomer',
            name='signal_sequence',
            field=cyano.history.HistoryForeignKey(related_name='protein_monomers', on_delete=django.db.models.deletion.SET_NULL, verbose_name='Sequence sequence', blank=True, to='cyano.SignalSequence', null=True),
        ),
        migrations.AddField(
            model_name='massspectrometryproteindetail',
            name='protein',
            field=cyano.history.HistoryForeignKey(related_name='protein_details', to='cyano.MassSpectrometryProtein'),
        ),
        migrations.AddField(
            model_name='disulfidebond',
            name='protein_monomer',
            field=cyano.history.HistoryForeignKey(related_name='disulfide_bonds', verbose_name='Protein monomer', to='cyano.ProteinMonomer'),
        ),
    ]
