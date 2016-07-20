# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models
import cyano.history
from django.conf import settings
import re
import django.core.validators
import django.db.models.deletion
import cyano.models
import datetime


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('contenttypes', '0002_remove_content_type_name'),
        ('auth', '0006_require_contenttypes_0002'),
    ]

    operations = [
        migrations.CreateModel(
            name='Basket',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('name', models.CharField(default='', verbose_name='Basket name', max_length=255)),
            ],
        ),
        migrations.CreateModel(
            name='BasketComponent',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('basket', cyano.history.HistoryForeignKey(verbose_name='In basket', related_name='components', to='cyano.Basket')),
            ],
        ),
        migrations.CreateModel(
            name='BindingSite',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('coordinate', models.PositiveIntegerField(verbose_name='Coordinate (nt)')),
                ('length', models.PositiveIntegerField(verbose_name='Length (nt)')),
                ('direction', models.CharField(choices=[('f', 'Forward'), ('r', 'Reverse')], verbose_name='Direction', max_length=10)),
            ],
            options={
                'verbose_name_plural': 'Binding sites',
                'verbose_name': 'Binding site',
                'ordering': ['coordinate', 'length'],
            },
        ),
        migrations.CreateModel(
            name='BiomassComposition',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('concentration', models.FloatField(verbose_name='Concentration (mmol gDCW<sup>-1</sup>)', validators=[django.core.validators.MinValueValidator(0)])),
            ],
            options={
                'verbose_name_plural': 'Biomass composition',
                'verbose_name': 'Biomass composition',
                'ordering': ['-concentration'],
            },
        ),
        migrations.CreateModel(
            name='Codon',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('sequence', models.CharField(verbose_name='Sequence', validators=[cyano.models.validate_dna_sequence, django.core.validators.MinLengthValidator(3)], max_length=3)),
            ],
            options={
                'verbose_name_plural': 'Codons',
                'verbose_name': 'Codon',
                'ordering': ['sequence'],
            },
        ),
        migrations.CreateModel(
            name='CoenzymeParticipant',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('coefficient', models.FloatField(validators=[django.core.validators.MinValueValidator(0)], null=True, verbose_name='Coefficient', blank=True)),
            ],
            options={
                'verbose_name_plural': 'Coenzyme participants',
                'verbose_name': 'Coenzyme participant',
                'ordering': [],
            },
        ),
        migrations.CreateModel(
            name='CrossReference',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('xid', models.CharField(verbose_name='External ID', max_length=255)),
                ('source', models.CharField(choices=[('ATCC', 'ATCC'), ('BiGG', 'BiGG'), ('BioCyc', 'BioCyc'), ('BioProject', 'BioProject'), ('CAS', 'CAS'), ('ChEBI', 'ChEBI'), ('CMR', 'CMR'), ('EC', 'EC'), ('GenBank', 'GenBank'), ('ISBN', 'ISBN'), ('KEGG', 'KEGG'), ('KNApSAcK', 'KNApSAcK'), ('LipidBank', 'LipidBank'), ('LIPIDMAPS', 'LIPIDMAPS'), ('PDB', 'PDB'), ('PDBCCD', 'PDBCCD'), ('PubChem', 'PubChem'), ('PubMed', 'PubMed'), ('RefSeq', 'RefSeq'), ('SABIO-RK', 'SABIO-RK'), ('SwissProt', 'SwissProt'), ('Taxonomy', 'Taxonomy'), ('ThreeDMET', 'ThreeDMET'), ('URL', 'URL')], verbose_name='Source', max_length=20)),
            ],
            options={
                'verbose_name_plural': 'Cross references',
                'verbose_name': 'Cross reference',
                'ordering': ['xid'],
            },
        ),
        migrations.CreateModel(
            name='DisulfideBond',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('residue_1', models.PositiveIntegerField(verbose_name='Residue-1')),
                ('residue_2', models.PositiveIntegerField(verbose_name='Residue-2')),
            ],
            options={
                'verbose_name_plural': 'Disulfide bonds',
                'verbose_name': 'Disulfide bond',
                'ordering': [],
            },
        ),
        migrations.CreateModel(
            name='DNAFootprint',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('length', models.PositiveIntegerField(null=True, verbose_name='Length (nt)', blank=True)),
                ('binding', models.CharField(choices=[('dsDNA', 'dsDNA'), ('ssDNA', 'ssDNA'), ('xsDNA', 'xsDNA')], verbose_name='Binding', blank=True, default='', max_length=10)),
                ('region', models.CharField(choices=[('dsDNA', 'dsDNA'), ('ssDNA', 'ssDNA'), ('xsDNA', 'xsDNA')], verbose_name='Region', blank=True, default='', max_length=10)),
            ],
            options={
                'verbose_name_plural': 'DNA footprints',
                'verbose_name': 'DNA footprint',
                'ordering': [],
            },
        ),
        migrations.CreateModel(
            name='Entry',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('wid', models.SlugField(verbose_name='WID', validators=[django.core.validators.RegexValidator(re.compile('^[-a-zA-Z0-9_]+\\Z', 32), "Enter a valid 'slug' consisting of letters, numbers, underscores or hyphens.", 'invalid')], max_length=150)),
                ('name', models.CharField(default='', verbose_name='Name', blank=True, max_length=255)),
                ('comments', models.TextField(default='', verbose_name='Comments', blank=True)),
            ],
            options={
                'verbose_name_plural': 'Entries',
                'listing': ['wid', 'name'],
                'permissions': (('view_normal', 'View entry'), ('view_delete', 'View deleted revisions of entry'), ('view_permission', 'View permissions of entry'), ('view_history', 'View older (not deleted) revisions of entry'), ('edit_permission', 'Allow modifying of permissions')),
                'verbose_name': 'Entry',
                'concrete_entry_model': False,
                'get_latest_by': 'createdDate',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'comments'],
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Comments', {'fields': ['comments']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
                'facet_fields': [],
                'ordering': ['wid'],
                'wid_unique': False,
            },
        ),
        migrations.CreateModel(
            name='EntryBasicTextData',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('value', models.TextField(default='', verbose_name='Value')),
            ],
            options={
                'verbose_name_plural': 'Entry basic text data',
                'verbose_name': 'Entry basic text data',
            },
        ),
        migrations.CreateModel(
            name='EntryBooleanData',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('value', models.BooleanField(verbose_name='Value')),
            ],
            options={
                'verbose_name_plural': 'Entry Boolean data',
                'verbose_name': 'Entry Boolean data',
                'ordering': ['value'],
            },
        ),
        migrations.CreateModel(
            name='EntryCharData',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('value', models.CharField(default='', verbose_name='Value', blank=True, max_length=255)),
                ('units', models.CharField(default='', verbose_name='Units', blank=True, max_length=255)),
            ],
            options={
                'verbose_name_plural': 'Entry char data',
                'verbose_name': 'Entry char data',
                'ordering': ['value', 'units'],
            },
        ),
        migrations.CreateModel(
            name='EntryFloatData',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('value', models.FloatField(verbose_name='Value')),
                ('units', models.CharField(default='', verbose_name='Units', blank=True, max_length=255)),
            ],
            options={
                'verbose_name_plural': 'Entry float data',
                'verbose_name': 'Entry float data',
                'ordering': ['value', 'units'],
            },
        ),
        migrations.CreateModel(
            name='EntryGroupObjectPermission',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='EntryPositiveFloatData',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('value', models.FloatField(verbose_name='Value', validators=[django.core.validators.MinValueValidator(0)])),
                ('units', models.CharField(default='', verbose_name='Units', blank=True, max_length=255)),
            ],
            options={
                'verbose_name_plural': 'Entry positive float data',
                'verbose_name': 'Entry positive float data',
                'ordering': ['value', 'units'],
            },
        ),
        migrations.CreateModel(
            name='EntryTextData',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('value', models.TextField(default='', verbose_name='Value', blank=True)),
                ('units', models.CharField(default='', verbose_name='Units', blank=True, max_length=255)),
            ],
            options={
                'verbose_name_plural': 'Entry text data',
                'verbose_name': 'Entry text data',
                'ordering': ['value', 'units'],
            },
        ),
        migrations.CreateModel(
            name='EntryUserObjectPermission',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='EnzymeParticipant',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
            ],
            options={
                'verbose_name_plural': 'Enzyme participants',
                'verbose_name': 'Enzyme participant',
                'ordering': [],
            },
        ),
        migrations.CreateModel(
            name='Evidence',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('value', models.TextField(default='', verbose_name='Value', blank=True)),
                ('units', models.CharField(null=True, verbose_name='Units', blank=True, max_length=255)),
                ('is_experimentally_constrained', models.BooleanField(verbose_name='Is experimentally <br/>constrained')),
                ('species', models.CharField(null=True, verbose_name='Species', blank=True, max_length=255)),
                ('media', models.CharField(null=True, verbose_name='Media', blank=True, max_length=255)),
                ('pH', models.FloatField(null=True, verbose_name='pH', blank=True)),
                ('temperature', models.FloatField(null=True, verbose_name='Temperature (C)', blank=True)),
                ('comments', models.TextField(default='', verbose_name='Comments', blank=True)),
            ],
            options={
                'verbose_name_plural': 'Evidence',
                'verbose_name': 'Evidence',
                'ordering': ['value', 'units'],
            },
        ),
        migrations.CreateModel(
            name='FeaturePosition',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('coordinate', models.PositiveIntegerField(verbose_name='Coordinate (nt)')),
                ('length', models.PositiveIntegerField(verbose_name='Length (nt)')),
                ('direction', models.CharField(choices=[('f', 'Forward'), ('r', 'Reverse')], verbose_name='Direction', max_length=10)),
            ],
            options={
                'verbose_name_plural': 'Feature Positions',
                'field_list': ['id', 'chromosome_feature', 'chromosome', 'coordinate', 'length'],
                'verbose_name': 'Feature Position',
                'fieldsets': [('Structure', {'fields': [{'verbose_name': 'Structure', 'name': 'structure'}, {'verbose_name': 'Structure Filter', 'name': 'structure_filter'}, {'verbose_name': 'Sequence', 'name': 'sequence'}, {'verbose_name': 'Genes', 'name': 'genes'}, {'verbose_name': 'Transcription units', 'name': 'transcription_units'}]})],
            },
        ),
        migrations.CreateModel(
            name='GlobalPermission',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
            ],
            options={
                'permissions': (('access_species', 'Can access any species'), ('create_mutant', 'Can create new species or mutants'), ('access_sbgn', 'Can access SBGN map')),
            },
        ),
        migrations.CreateModel(
            name='GroupProfile',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('description', models.CharField(default='', verbose_name='Group description', blank=True, max_length=255)),
                ('group', models.OneToOneField(related_name='profile', to='auth.Group')),
            ],
        ),
        migrations.CreateModel(
            name='Homolog',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('xid', models.CharField(verbose_name='External ID', max_length=255)),
                ('species', models.CharField(choices=[('B. subtilis', 'B. subtilis'), ('E. coli', 'E. coli'), ('M. hyopneumoniae', 'M. hyopneumoniae'), ('M. mobile', 'M. mobile'), ('M. pneumoniae', 'M. pneumoniae'), ('S. coelicolor', 'S. coelicolor'), ('S. oneidensis', 'S. oneidensis')], verbose_name='Species', max_length=20)),
                ('evidence', cyano.history.HistoryManyToManyField(verbose_name='Evidence', blank=True, to='cyano.Evidence')),
            ],
            options={
                'verbose_name_plural': 'Homologs',
                'verbose_name': 'Homolog',
                'ordering': ['xid'],
            },
        ),
        migrations.CreateModel(
            name='Kinetics',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('rate_law', models.CharField(default='', verbose_name='Rate law', blank=True, max_length=255)),
                ('km', models.CharField(validators=[django.core.validators.RegexValidator('^([0-9\\.]+)(, [0-9\\.]+)*$')], verbose_name='K<sub>m</sub> (&mu;M)', blank=True, max_length=255)),
                ('vmax', models.FloatField(validators=[django.core.validators.MinValueValidator(0)], null=True, verbose_name='V<sub>max</sub>', blank=True)),
                ('vmax_unit', models.CharField(choices=[('1/min', '1/min'), ('U/mg', 'U/mg')], verbose_name='V<sub>max</sub> Unit', blank=True, max_length=255)),
                ('evidence', cyano.history.HistoryManyToManyField(verbose_name='Evidence', blank=True, to='cyano.Evidence')),
            ],
            options={
                'verbose_name_plural': 'Kinetics',
                'verbose_name': 'Kinetics',
                'ordering': [],
            },
        ),
        migrations.CreateModel(
            name='MassSpectrometryProteinDetail',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
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
                'verbose_name_plural': 'Mass Spectrometry Protein Details',
                'field_list': ['id', 'chromosome_feature', 'chromosome', 'coordinate', 'length'],
                'verbose_name': 'Mass Spectrometry Protein Detail',
                'fieldsets': [('Mass Spectrometry Details', {'fields': ['sequence', 'sequence_ptm', 'proteotypic', 'zscore', 'delta_mass', 'mass', 'charge', 'retention_time', 'theoretical_mass', 'missed_cleavages']})],
            },
        ),
        migrations.CreateModel(
            name='MediaComposition',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('concentration', models.FloatField(verbose_name='Concentration (mM)', validators=[django.core.validators.MinValueValidator(0)])),
                ('is_diffused', models.BooleanField(verbose_name='Is diffused')),
                ('evidence', cyano.history.HistoryManyToManyField(verbose_name='Evidence', blank=True, to='cyano.Evidence')),
            ],
            options={
                'verbose_name_plural': 'Media composition',
                'verbose_name': 'Media composition',
                'ordering': ['-concentration'],
            },
        ),
        migrations.CreateModel(
            name='MetaboliteMapCoordinate',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('x', models.FloatField(verbose_name='X')),
                ('y', models.FloatField(verbose_name='Y')),
            ],
            options={
                'verbose_name_plural': 'Metabolite map coordinates',
                'verbose_name': 'Metabolite map coordinate',
                'ordering': ['x', 'y', 'compartment'],
            },
        ),
        migrations.CreateModel(
            name='ModificationReaction',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('position', models.PositiveIntegerField(null=True, verbose_name='Position', blank=True)),
                ('evidence', cyano.history.HistoryManyToManyField(verbose_name='Evidence', blank=True, to='cyano.Evidence')),
            ],
            options={
                'verbose_name_plural': 'Protein monomers',
                'verbose_name': 'Protein monomer',
                'ordering': [],
            },
        ),
        migrations.CreateModel(
            name='ProstheticGroupParticipant',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('coefficient', models.PositiveIntegerField(null=True, verbose_name='Coefficient', blank=True)),
                ('evidence', cyano.history.HistoryManyToManyField(verbose_name='Evidence', blank=True, to='cyano.Evidence')),
            ],
            options={
                'verbose_name_plural': 'Prosthetic group participants',
                'verbose_name': 'Prosthetic group participant',
                'ordering': [],
            },
        ),
        migrations.CreateModel(
            name='ProteinComparison',
            fields=[
                ('wid', models.AutoField(primary_key=True, serialize=False)),
                ('protein_a', models.IntegerField()),
                ('protein_b', models.IntegerField()),
                ('equal_score', models.FloatField()),
            ],
        ),
        migrations.CreateModel(
            name='ProteinComplexBiosythesisParticipant',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('residue', models.PositiveIntegerField(null=True, verbose_name='Residue', blank=True)),
                ('coefficient', models.FloatField(verbose_name='Coefficient')),
                ('evidence', cyano.history.HistoryManyToManyField(verbose_name='Evidence', blank=True, to='cyano.Evidence')),
            ],
            options={
                'verbose_name_plural': 'Protein complex biosythesis participants',
                'verbose_name': 'Protein complex biosythesis participant',
                'ordering': [],
            },
        ),
        migrations.CreateModel(
            name='ReactionMapCoordinate',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('path', models.TextField(verbose_name='Path')),
                ('value_x', models.FloatField(verbose_name='Label-X')),
                ('value_y', models.FloatField(verbose_name='Label-Y')),
                ('label_x', models.FloatField(verbose_name='Value-X')),
                ('label_y', models.FloatField(verbose_name='Value-Y')),
            ],
            options={
                'verbose_name_plural': 'Reaction map coordinates',
                'verbose_name': 'Reaction map coordinate',
                'ordering': ['value_x', 'value_y', 'label_x', 'label_y'],
            },
        ),
        migrations.CreateModel(
            name='ReactionStoichiometryParticipant',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('coefficient', models.FloatField(verbose_name='Coefficient')),
                ('evidence', cyano.history.HistoryManyToManyField(verbose_name='Evidence', blank=True, to='cyano.Evidence')),
            ],
            options={
                'verbose_name_plural': 'Molecule coefficient compartments',
                'verbose_name': 'Molecule coefficient compartment',
                'ordering': [],
            },
        ),
        migrations.CreateModel(
            name='Revision',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('object_id', models.IntegerField(db_index=True, verbose_name='Current version primary key')),
                ('action', models.CharField(choices=[('I', 'Insert'), ('U', 'Update'), ('D', 'Delete'), ('X', 'Unknown')], max_length=1)),
                ('new_data', models.TextField(blank=True)),
                ('content_type', cyano.history.HistoryForeignKey(to='contenttypes.ContentType')),
            ],
            options={
                'ordering': ['detail__date'],
            },
        ),
        migrations.CreateModel(
            name='RevisionDetail',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('date', models.DateTimeField(default=datetime.datetime.now, verbose_name='Modification date')),
                ('reason', models.TextField(default='', verbose_name='Reason for edit', blank=True)),
            ],
        ),
        migrations.CreateModel(
            name='SignalSequence',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('type', models.CharField(choices=[('lipoprotein', 'Lipoprotein'), ('secretory', 'Secretory')], verbose_name='Type', max_length=20)),
                ('location', models.CharField(choices=[('N', 'N-terminus'), ('C', 'C-terminus')], verbose_name='Location', max_length=1)),
                ('length', models.PositiveIntegerField(verbose_name='Length (nt)')),
                ('evidence', cyano.history.HistoryManyToManyField(verbose_name='Evidence', blank=True, to='cyano.Evidence')),
            ],
            options={
                'verbose_name_plural': 'Signal sequences',
                'verbose_name': 'Signal sequence',
                'ordering': ['type', 'location', 'length'],
            },
        ),
        migrations.CreateModel(
            name='Synonym',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('name', models.CharField(verbose_name='Name', max_length=255)),
            ],
            options={
                'verbose_name_plural': 'Synonyms',
                'verbose_name': 'Synonym',
                'ordering': ['name'],
            },
        ),
        migrations.CreateModel(
            name='TableMeta',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('table_name', models.CharField(unique=True, verbose_name='Name of the table', max_length=255)),
                ('model_name', models.CharField(unique=True, verbose_name='Name of the model associated with the table', max_length=255)),
            ],
        ),
        migrations.CreateModel(
            name='UserProfile',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('affiliation', models.CharField(default='', verbose_name='Affiliation', blank=True, max_length=255)),
                ('website', models.URLField(default='', verbose_name='Website', blank=True, max_length=255)),
                ('phone', models.CharField(default='', verbose_name='Phone', blank=True, max_length=255)),
                ('address', models.CharField(default='', verbose_name='Address', blank=True, max_length=255)),
                ('city', models.CharField(default='', verbose_name='City', blank=True, max_length=255)),
                ('state', models.CharField(default='', verbose_name='State', blank=True, max_length=255)),
                ('zip', models.CharField(default='', verbose_name='Zip', blank=True, max_length=255)),
                ('country', models.CharField(default='', verbose_name='Country', blank=True, max_length=255)),
                ('force_password_change', models.BooleanField(default=False)),
                ('user', models.OneToOneField(related_name='profile', to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'verbose_name_plural': 'User profiles',
                'verbose_name': 'User profile',
                'get_latest_by': 'user__date_joined',
                'ordering': ['user__last_name', 'user__first_name'],
            },
        ),
        migrations.CreateModel(
            name='Species',
            fields=[
                ('parent_ptr_entry', models.OneToOneField(to='cyano.Entry', verbose_name='Entry', serialize=False, related_name='child_ptr_species', parent_link=True, primary_key=True)),
                ('genetic_code', models.CharField(choices=[('1', 'Standard'), ('2', 'Vertebrate'), ('3', 'Yeast'), ('4', 'Mold, protozoa, coelenterate mitochondria, mycoplasma, and spiroplasma'), ('5', 'Invertebrate mitochondria'), ('6', 'Ciliate, dasycladacean and hexamita'), ('9', 'Echinoderm and flatworm mitochondria'), ('10', 'Euplotid'), ('11', 'Bacteria, archaea and plant plastids'), ('12', 'Alternative yeast'), ('13', 'Ascidian mitochondria'), ('14', 'Alternative flatworm mitochondria'), ('15', 'Blepharisma'), ('16', 'Chlorophycean mitochondria'), ('21', 'Trematode mitochondria'), ('22', 'Scenedesmus obliquus mitochondria'), ('23', 'Thraustochytrium mitochondria'), ('24', 'Pterobranchia mitochondria')], verbose_name='Genetic code', max_length=50)),
            ],
            options={
                'verbose_name_plural': 'Species',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'publication_references', 'cross_references', 'genetic_code', 'comments'],
                'verbose_name': 'Species',
                'concrete_entry_model': True,
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Physiology', {'fields': ['genetic_code']}), ('Comments', {'fields': ['comments']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
                'facet_fields': [],
                'wid_unique': True,
            },
            bases=('cyano.entry',),
        ),
        migrations.CreateModel(
            name='SpeciesComponent',
            fields=[
                ('entry_ptr', models.OneToOneField(primary_key=True, serialize=False, parent_link=True, to='cyano.Entry', auto_created=True)),
            ],
            options={
                'verbose_name_plural': 'Species components',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'comments', 'publication_references'],
                'verbose_name': 'Species component',
                'concrete_entry_model': False,
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
                'facet_fields': ['type'],
                'wid_unique': False,
            },
            bases=('cyano.entry',),
        ),
        migrations.AddField(
            model_name='revisiondetail',
            name='user',
            field=cyano.history.HistoryForeignKey(editable=False, verbose_name='Modified by', related_name='+', to='cyano.UserProfile'),
        ),
        migrations.AddField(
            model_name='revision',
            name='detail',
            field=cyano.history.HistoryForeignKey(editable=False, verbose_name='Details about this revision', related_name='revisions', to='cyano.RevisionDetail'),
        ),
        migrations.AddField(
            model_name='enzymeparticipant',
            name='evidence',
            field=cyano.history.HistoryManyToManyField(verbose_name='Evidence', blank=True, to='cyano.Evidence'),
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
            field=cyano.history.HistoryManyToManyField(verbose_name='Evidence', blank=True, to='cyano.Evidence'),
        ),
        migrations.AddField(
            model_name='entrypositivefloatdata',
            name='evidence',
            field=cyano.history.HistoryManyToManyField(verbose_name='Evidence', blank=True, to='cyano.Evidence'),
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
            field=cyano.history.HistoryManyToManyField(verbose_name='Evidence', blank=True, to='cyano.Evidence'),
        ),
        migrations.AddField(
            model_name='entrychardata',
            name='evidence',
            field=cyano.history.HistoryManyToManyField(verbose_name='Evidence', blank=True, to='cyano.Evidence'),
        ),
        migrations.AddField(
            model_name='entrybooleandata',
            name='evidence',
            field=cyano.history.HistoryManyToManyField(verbose_name='Evidence', blank=True, to='cyano.Evidence'),
        ),
        migrations.AddField(
            model_name='entry',
            name='model_type',
            field=cyano.history.HistoryForeignKey(to='cyano.TableMeta'),
        ),
        migrations.AddField(
            model_name='entry',
            name='synonyms',
            field=cyano.history.HistoryManyToManyField(verbose_name='Synonyms', blank=True, related_name='entry', to='cyano.Synonym'),
        ),
        migrations.AddField(
            model_name='dnafootprint',
            name='evidence',
            field=cyano.history.HistoryManyToManyField(verbose_name='Evidence', blank=True, to='cyano.Evidence'),
        ),
        migrations.AddField(
            model_name='disulfidebond',
            name='evidence',
            field=cyano.history.HistoryManyToManyField(verbose_name='Evidence', blank=True, to='cyano.Evidence'),
        ),
        migrations.AlterUniqueTogether(
            name='crossreference',
            unique_together=set([('xid', 'source')]),
        ),
        migrations.AddField(
            model_name='coenzymeparticipant',
            name='evidence',
            field=cyano.history.HistoryManyToManyField(verbose_name='Evidence', blank=True, to='cyano.Evidence'),
        ),
        migrations.AddField(
            model_name='codon',
            name='evidence',
            field=cyano.history.HistoryManyToManyField(verbose_name='Evidence', blank=True, to='cyano.Evidence'),
        ),
        migrations.AddField(
            model_name='biomasscomposition',
            name='evidence',
            field=cyano.history.HistoryManyToManyField(verbose_name='Evidence', blank=True, to='cyano.Evidence'),
        ),
        migrations.AddField(
            model_name='bindingsite',
            name='evidence',
            field=cyano.history.HistoryManyToManyField(verbose_name='Evidence', blank=True, to='cyano.Evidence'),
        ),
        migrations.AddField(
            model_name='basket',
            name='user',
            field=cyano.history.HistoryForeignKey(verbose_name='Users baskets', related_name='baskets', to='cyano.UserProfile'),
        ),
        migrations.CreateModel(
            name='ChromosomeFeature',
            fields=[
                ('parent_ptr_species_component', models.OneToOneField(to='cyano.SpeciesComponent', verbose_name='Species component', serialize=False, related_name='child_ptr_chromosome_feature', parent_link=True, primary_key=True)),
            ],
            options={
                'verbose_name_plural': 'Chromosome features',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'comments', 'publication_references'],
                'verbose_name': 'Chromosome feature',
                'concrete_entry_model': True,
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), {'inline': 'positions'}, ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
                'facet_fields': ['type'],
                'wid_unique': False,
            },
            bases=('cyano.speciescomponent',),
        ),
        migrations.CreateModel(
            name='Compartment',
            fields=[
                ('parent_ptr_species_component', models.OneToOneField(to='cyano.SpeciesComponent', verbose_name='Species component', serialize=False, related_name='child_ptr_compartment', parent_link=True, primary_key=True)),
            ],
            options={
                'verbose_name_plural': 'Compartments',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'comments', 'publication_references'],
                'verbose_name': 'Compartment',
                'concrete_entry_model': True,
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Content', {'fields': [{'verbose_name': 'Metabolites (mM)', 'name': 'biomass_compositions'}, {'verbose_name': 'Protein monomers', 'name': 'protein_monomers'}, {'verbose_name': 'Protein complexes', 'name': 'protein_complexes'}]}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
                'facet_fields': ['type'],
                'wid_unique': False,
            },
            bases=('cyano.speciescomponent',),
        ),
        migrations.CreateModel(
            name='MassSpectrometryJob',
            fields=[
                ('parent_ptr_species_component', models.OneToOneField(to='cyano.SpeciesComponent', verbose_name='Species component', serialize=False, related_name='child_ptr_mass_spectrometry', parent_link=True, primary_key=True)),
            ],
            options={
                'verbose_name_plural': 'Mass Spectrometry Jobs',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'comments', 'publication_references'],
                'verbose_name': 'Mass Spectrometry Job',
                'concrete_entry_model': True,
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Mass Spectrometry', {'fields': [{'verbose_name': 'Target Peptides', 'name': 'target_peptide'}, {'verbose_name': 'Decoy Peptides', 'name': 'decoy_peptide'}, {'verbose_name': 'Related Proteins', 'name': 'related_proteins'}]}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
                'facet_fields': ['type'],
                'wid_unique': False,
            },
            bases=('cyano.speciescomponent',),
        ),
        migrations.CreateModel(
            name='Molecule',
            fields=[
                ('parent_ptr_species_component', models.OneToOneField(to='cyano.SpeciesComponent', verbose_name='Species component', serialize=False, related_name='child_ptr_molecule', parent_link=True, primary_key=True)),
            ],
            options={
                'verbose_name_plural': 'Molecules',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'comments', 'publication_references'],
                'verbose_name': 'Molecule',
                'concrete_entry_model': False,
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Function', {'fields': [{'verbose_name': 'Reaction participant', 'name': 'reaction_stoichiometry_participants'}, {'verbose_name': 'Complex subunit', 'name': 'protein_complex_biosythesis_participants'}]}), ('Parameters', {'fields': ['parameters']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
                'facet_fields': ['type'],
                'wid_unique': False,
            },
            bases=('cyano.speciescomponent',),
        ),
        migrations.CreateModel(
            name='Note',
            fields=[
                ('parent_ptr_species_component', models.OneToOneField(to='cyano.SpeciesComponent', verbose_name='Species component', serialize=False, related_name='child_ptr_note', parent_link=True, primary_key=True)),
            ],
            options={
                'verbose_name_plural': 'Notes',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'comments', 'publication_references'],
                'verbose_name': 'Note',
                'concrete_entry_model': True,
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
                'facet_fields': ['type'],
                'wid_unique': False,
            },
            bases=('cyano.speciescomponent',),
        ),
        migrations.CreateModel(
            name='Parameter',
            fields=[
                ('parent_ptr_species_component', models.OneToOneField(to='cyano.SpeciesComponent', verbose_name='Species component', serialize=False, related_name='child_ptr_parameter', parent_link=True, primary_key=True)),
            ],
            options={
                'verbose_name_plural': 'Misc. parameters',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'value', 'reactions', 'molecules', 'state', 'process', 'comments', 'publication_references'],
                'verbose_name': 'Misc. parameter',
                'concrete_entry_model': True,
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Value', {'fields': ['value']}), ('Associations', {'fields': ['reactions', 'molecules', 'state', 'process']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
                'facet_fields': ['type', 'reactions', 'molecules', 'state', 'process'],
                'wid_unique': False,
            },
            bases=('cyano.speciescomponent',),
        ),
        migrations.CreateModel(
            name='Pathway',
            fields=[
                ('parent_ptr_species_component', models.OneToOneField(to='cyano.SpeciesComponent', verbose_name='Species component', serialize=False, related_name='child_ptr_pathway', parent_link=True, primary_key=True)),
            ],
            options={
                'verbose_name_plural': 'Pathways',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'comments', 'publication_references', 'type'],
                'verbose_name': 'Pathway',
                'concrete_entry_model': True,
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Navigator', {'fields': [{'verbose_name': 'Navigate to', 'name': 'navigator'}]}), ('Reactions', {'fields': [{'verbose_name': 'Reactions', 'name': 'reaction_map'}, 'reactions']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
                'facet_fields': ['type'],
                'group_field': 'type',
                'wid_unique': True,
            },
            bases=('cyano.speciescomponent',),
        ),
        migrations.CreateModel(
            name='Process',
            fields=[
                ('parent_ptr_species_component', models.OneToOneField(to='cyano.SpeciesComponent', verbose_name='Species component', serialize=False, related_name='child_ptr_process', parent_link=True, primary_key=True)),
                ('initialization_order', models.PositiveIntegerField(verbose_name='Initialization order')),
            ],
            options={
                'verbose_name_plural': 'Processes',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'initialization_order', 'comments', 'publication_references'],
                'verbose_name': 'Process',
                'concrete_entry_model': True,
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Implementation', {'fields': ['initialization_order']}), ('Reactions', {'fields': [{'verbose_name': 'Chemical reactions', 'name': 'reactions'}, {'verbose_name': 'Complex formation reactions', 'name': 'formed_complexes'}]}), ('Parameters', {'fields': ['parameters']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
                'facet_fields': ['type'],
                'wid_unique': False,
            },
            bases=('cyano.speciescomponent',),
        ),
        migrations.CreateModel(
            name='PublicationReference',
            fields=[
                ('parent_ptr_species_component', models.OneToOneField(to='cyano.SpeciesComponent', verbose_name='Species component', serialize=False, related_name='child_ptr_publicationreference', parent_link=True, primary_key=True)),
                ('authors', models.TextField(default='', verbose_name='Author(s)', blank=True)),
                ('editors', models.TextField(default='', verbose_name='Editor(s)', blank=True)),
                ('year', models.PositiveIntegerField(null=True, verbose_name='Year', blank=True)),
                ('title', models.TextField(default='', verbose_name='Title', blank=True)),
                ('publication', models.CharField(default='', verbose_name='Publication', blank=True, max_length=255)),
                ('publisher', models.CharField(default='', verbose_name='Publisher', blank=True, max_length=255)),
                ('volume', models.CharField(default='', verbose_name='Volume', blank=True, max_length=255)),
                ('issue', models.CharField(default='', verbose_name='Issue', blank=True, max_length=255)),
                ('pages', models.CharField(default='', verbose_name='Page(s)', blank=True, max_length=255)),
            ],
            options={
                'verbose_name_plural': 'Publication References',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'authors', 'editors', 'year', 'title', 'publication', 'publisher', 'volume', 'issue', 'pages', 'comments', 'publication_references'],
                'verbose_name': 'Publication Reference',
                'concrete_entry_model': True,
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Citation', {'fields': [{'verbose_name': 'Citation', 'name': 'citation'}]}), ('Cited by', {'fields': [{'verbose_name': 'Cited by', 'name': 'referenced_entries'}]}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
                'facet_fields': ['type', 'year', 'publication'],
                'wid_unique': True,
            },
            bases=('cyano.speciescomponent',),
        ),
        migrations.CreateModel(
            name='Reaction',
            fields=[
                ('parent_ptr_species_component', models.OneToOneField(to='cyano.SpeciesComponent', verbose_name='Species component', serialize=False, related_name='child_ptr_reaction', parent_link=True, primary_key=True)),
                ('direction', models.CharField(choices=[('f', 'Forward'), ('b', 'Backward'), ('r', 'Reversible')], verbose_name='Direction', max_length=1)),
                ('is_spontaneous', models.BooleanField(verbose_name='Is spontaneous (pH 7.5, 25C, <i>I</i> = 0)')),
                ('delta_g', models.FloatField(null=True, verbose_name='&Delta;G (pH 7.5, 25C, <i>I</i> = 0; kJ mol<sup>-1</sup>)', blank=True)),
                ('coenzymes', cyano.history.HistoryManyToManyField(verbose_name='Coenzymes', blank=True, related_name='reactions', to='cyano.CoenzymeParticipant')),
                ('enzyme', cyano.history.HistoryForeignKey(null=True, verbose_name='Enzyme', blank=True, related_name='reactions', to='cyano.EnzymeParticipant', on_delete=django.db.models.deletion.SET_NULL)),
                ('keq', cyano.history.HistoryForeignKey(null=True, verbose_name='K<sub>eq</sub>', blank=True, related_name='+', to='cyano.EntryPositiveFloatData', on_delete=django.db.models.deletion.SET_NULL)),
                ('kinetics_backward', cyano.history.HistoryForeignKey(null=True, verbose_name='Backward kinetics', blank=True, related_name='reactions_backward', to='cyano.Kinetics', on_delete=django.db.models.deletion.SET_NULL)),
                ('kinetics_forward', cyano.history.HistoryForeignKey(null=True, verbose_name='Forward kinetics', blank=True, related_name='reactions_forward', to='cyano.Kinetics', on_delete=django.db.models.deletion.SET_NULL)),
                ('map_coordinates', cyano.history.HistoryManyToManyField(verbose_name='Map coordinates', blank=True, related_name='reactions', to='cyano.ReactionMapCoordinate')),
                ('modification', cyano.history.HistoryForeignKey(null=True, verbose_name='Modification', blank=True, related_name='reactions', to='cyano.ModificationReaction', on_delete=django.db.models.deletion.SET_NULL)),
                ('optimal_ph', cyano.history.HistoryForeignKey(null=True, verbose_name='Optimal pH', blank=True, related_name='+', to='cyano.EntryPositiveFloatData', on_delete=django.db.models.deletion.SET_NULL)),
                ('optimal_temperature', cyano.history.HistoryForeignKey(null=True, verbose_name='Optimal temperature', blank=True, related_name='+', to='cyano.EntryFloatData', on_delete=django.db.models.deletion.SET_NULL)),
                ('pathways', cyano.history.HistoryManyToManyField(verbose_name='Pathways', blank=True, related_name='reactions', to='cyano.Pathway')),
                ('processes', cyano.history.HistoryForeignKey(null=True, verbose_name='Process', blank=True, related_name='reactions', to='cyano.Process', on_delete=django.db.models.deletion.SET_NULL)),
            ],
            options={
                'verbose_name_plural': 'Reactions',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'stoichiometry', 'direction', 'modification', 'enzyme', 'coenzymes', 'optimal_ph', 'optimal_temperature', 'is_spontaneous', 'delta_g', 'keq', 'kinetics_forward', 'kinetics_backward', 'pathways', 'processes', 'states', 'map_coordinates', 'comments', 'publication_references'],
                'verbose_name': 'Reaction',
                'concrete_entry_model': True,
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Reaction', {'fields': ['stoichiometry', 'modification']}), ('Catalysis', {'fields': ['enzyme', 'coenzymes', 'optimal_ph', 'optimal_temperature']}), ('Energetics', {'fields': ['is_spontaneous', 'delta_g', 'keq']}), ('Kinetics', {'fields': ['kinetics_forward', 'kinetics_backward']}), ('Parameters', {'fields': ['parameters']}), ('Associations', {'fields': ['pathways', 'processes', 'states']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
                'facet_fields': ['type', 'direction', 'enzyme__protein', 'coenzymes__metabolite', 'is_spontaneous', 'pathways', 'processes', 'states'],
                'wid_unique': False,
            },
            bases=('cyano.speciescomponent',),
        ),
        migrations.CreateModel(
            name='State',
            fields=[
                ('parent_ptr_species_component', models.OneToOneField(to='cyano.SpeciesComponent', verbose_name='Species component', serialize=False, related_name='child_ptr_state', parent_link=True, primary_key=True)),
            ],
            options={
                'verbose_name_plural': 'States',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'comments', 'publication_references'],
                'verbose_name': 'State',
                'concrete_entry_model': True,
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Reactions', {'fields': ['reactions']}), ('Parameters', {'fields': ['parameters']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
                'facet_fields': ['type'],
                'wid_unique': False,
            },
            bases=('cyano.speciescomponent',),
        ),
        migrations.CreateModel(
            name='TranscriptionalRegulation',
            fields=[
                ('parent_ptr_species_component', models.OneToOneField(to='cyano.SpeciesComponent', verbose_name='Species component', serialize=False, related_name='child_ptr_transcriptional_regulation', parent_link=True, primary_key=True)),
                ('activity', cyano.history.HistoryForeignKey(verbose_name='Fold-change activity', related_name='+', to='cyano.EntryPositiveFloatData')),
                ('affinity', cyano.history.HistoryForeignKey(null=True, verbose_name='Affinity', blank=True, related_name='+', to='cyano.EntryPositiveFloatData', on_delete=django.db.models.deletion.SET_NULL)),
                ('binding_site', cyano.history.HistoryForeignKey(null=True, verbose_name='Binding site', blank=True, related_name='transcriptional_regulations', to='cyano.BindingSite', on_delete=django.db.models.deletion.SET_NULL)),
            ],
            options={
                'verbose_name_plural': 'Transcriptional regulation',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'transcription_unit', 'transcription_factor', 'binding_site', 'affinity', 'activity', 'comments', 'publication_references'],
                'verbose_name': 'Transcriptional regulation',
                'concrete_entry_model': True,
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Regulation', {'fields': ['transcription_unit', 'transcription_factor', 'binding_site', 'affinity', 'activity']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
                'facet_fields': ['type', 'transcription_unit', 'transcription_factor'],
                'wid_unique': False,
            },
            bases=('cyano.speciescomponent',),
        ),
        migrations.CreateModel(
            name='Type',
            fields=[
                ('parent_ptr_species_component', models.OneToOneField(to='cyano.SpeciesComponent', verbose_name='Species component', serialize=False, related_name='child_ptr_type', parent_link=True, primary_key=True)),
            ],
            options={
                'verbose_name_plural': 'Types',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'parent', 'comments', 'publication_references'],
                'verbose_name': 'Type',
                'concrete_entry_model': True,
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type', 'parent', 'children', 'members']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
                'facet_fields': ['type', 'parent'],
                'wid_unique': True,
            },
            bases=('cyano.speciescomponent',),
        ),
        migrations.AddField(
            model_name='speciescomponent',
            name='cross_references',
            field=cyano.history.HistoryManyToManyField(verbose_name='Cross references', blank=True, related_name='cross_referenced_components', to='cyano.CrossReference'),
        ),
        migrations.AddField(
            model_name='speciescomponent',
            name='parent',
            field=cyano.history.HistoryForeignKey(null=True, verbose_name='Parent', blank=True, related_name='children', to='cyano.SpeciesComponent', on_delete=django.db.models.deletion.SET_NULL),
        ),
        migrations.AddField(
            model_name='speciescomponent',
            name='species',
            field=cyano.history.HistoryForeignKey(verbose_name='Species', related_name='components', to='cyano.Species'),
        ),
        migrations.AddField(
            model_name='species',
            name='cross_references',
            field=cyano.history.HistoryManyToManyField(verbose_name='Cross references', blank=True, related_name='cross_referenced_species', to='cyano.CrossReference'),
        ),
        migrations.AddField(
            model_name='evidence',
            name='species_component',
            field=cyano.history.HistoryManyToManyField(verbose_name='Species component', blank=True, related_name='_evidence_species_component_+', to='cyano.SpeciesComponent'),
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
            field=cyano.history.HistoryForeignKey(verbose_name='component', related_name='+', to='cyano.SpeciesComponent'),
        ),
        migrations.AddField(
            model_name='basketcomponent',
            name='species',
            field=cyano.history.HistoryForeignKey(verbose_name='Species component belongs to', related_name='+', to='cyano.Species'),
        ),
        migrations.CreateModel(
            name='Gene',
            fields=[
                ('parent_ptr_molecule', models.OneToOneField(to='cyano.Molecule', verbose_name='Molecule', serialize=False, related_name='child_ptr_gene', parent_link=True, primary_key=True)),
                ('symbol', models.CharField(default='', verbose_name='Symbol', blank=True, max_length=255)),
                ('coordinate', models.PositiveIntegerField(verbose_name='Coordinate (nt)')),
                ('length', models.PositiveIntegerField(verbose_name='Length (nt)')),
                ('direction', models.CharField(choices=[('f', 'Forward'), ('r', 'Reverse')], verbose_name='Direction', max_length=10)),
            ],
            options={
                'verbose_name_plural': 'Genes',
                'listing': ['wid', 'symbol'],
                'verbose_name': 'Gene',
                'concrete_entry_model': True,
                'field_list': ['id', 'wid', 'name', 'symbol', 'synonyms', 'cross_references', 'homologs', 'type', 'chromosome', 'coordinate', 'length', 'direction', 'is_essential', 'expression', 'half_life', 'codons', 'amino_acid', 'comments', 'publication_references'],
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'symbol', 'synonyms', {'verbose_name': 'Protein product', 'name': 'protein_monomers'}, 'cross_references', 'homologs']}), ('Classification', {'fields': ['type']}), ('Structure', {'fields': [{'verbose_name': 'Structure', 'name': 'structure'}, {'verbose_name': 'Structure Filter', 'name': 'structure_filter'}, {'verbose_name': 'Sequence', 'name': 'sequence'}, {'verbose_name': 'Transcription unit', 'name': 'transcription_units'}, {'verbose_name': 'Empirical formula (pH 7.5)', 'name': 'empirical_formula'}, {'verbose_name': 'Molecular weight (pH 7.5; Da)', 'name': 'molecular_weight'}]}), ('Functional genomics', {'fields': ['is_essential', 'expression', 'half_life', 'codons', 'amino_acid', {'verbose_name': 'Extinction coefficient <br/>(260 nm, 25C, pH 7.0)', 'name': 'extinction_coefficient'}, {'verbose_name': 'pI', 'name': 'pi'}]}), ('Function', {'fields': [{'verbose_name': 'Reaction participant', 'name': 'reaction_stoichiometry_participants'}, {'verbose_name': 'Complex subunit', 'name': 'protein_complex_biosythesis_participants'}]}), ('Parameters', {'fields': ['parameters']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
                'facet_fields': ['type', 'chromosome', 'direction', 'is_essential', 'amino_acid'],
                'wid_unique': False,
            },
            bases=('cyano.molecule',),
        ),
        migrations.CreateModel(
            name='Genome',
            fields=[
                ('parent_ptr_molecule', models.OneToOneField(to='cyano.Molecule', verbose_name='Molecule', serialize=False, related_name='child_ptr_genome', parent_link=True, primary_key=True)),
                ('sequence', models.TextField(validators=[cyano.models.validate_dna_sequence], default='', verbose_name='Sequence', blank=True)),
                ('length', models.PositiveIntegerField(verbose_name='Length (nt)')),
            ],
            options={
                'verbose_name_plural': 'Genome',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'sequence', 'length', 'comments', 'publication_references'],
                'verbose_name': 'Genome',
                'concrete_entry_model': True,
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Sequence', {'fields': [{'verbose_name': 'Structure Filter', 'name': 'structure_filter'}, {'verbose_name': 'Structure', 'name': 'structure'}, {'verbose_name': 'Extinction coefficient <br/>(260 nm, 25C, pH 7.0)', 'name': 'extinction_coefficient'}, {'verbose_name': 'pI', 'name': 'pi'}]}), ('Features', {'fields': ['genes', {'verbose_name': 'Other features', 'name': 'features'}]}), ('Function', {'fields': [{'verbose_name': 'Reaction participant', 'name': 'reaction_stoichiometry_participants'}, {'verbose_name': 'Complex subunit', 'name': 'protein_complex_biosythesis_participants'}]}), ('Parameters', {'fields': ['parameters']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
                'facet_fields': ['type'],
                'wid_unique': False,
            },
            bases=('cyano.molecule',),
        ),
        migrations.CreateModel(
            name='Metabolite',
            fields=[
                ('parent_ptr_molecule', models.OneToOneField(to='cyano.Molecule', verbose_name='Molecule', serialize=False, related_name='child_ptr_metabolite', parent_link=True, primary_key=True)),
                ('traditional_name', models.TextField(default='', verbose_name='Traditional name', blank=True)),
                ('iupac_name', models.TextField(default='', verbose_name='IUPAC name', blank=True)),
                ('empirical_formula', models.TextField(verbose_name='Empirical formula (pH 7.5)', validators=[django.core.validators.RegexValidator(message='Invalid empirical formula', regex='^([A-Z][a-z]*[0-9]*)+$')])),
                ('smiles', models.TextField(default='', verbose_name='SMILES (pH 7.5)', blank=True)),
                ('charge', models.IntegerField(verbose_name='Charge (pH 7.5)')),
                ('is_hydrophobic', models.BooleanField(verbose_name='Is hydrophobic')),
                ('volume', models.FloatField(validators=[django.core.validators.MinValueValidator(0)], null=True, verbose_name='van der Waals volume <br/>(pH 7.5; &#8491;<sup>3</sup> molecule<sup>-1</sup>)', blank=True)),
                ('deltag_formation', models.FloatField(null=True, verbose_name="&Delta;<sub>f</sub>G<sup>'o</sup> (pH 7.5, 25C, I = 0; kJ mol<sup>-1</sup>)", blank=True)),
                ('pka', models.FloatField(validators=[django.core.validators.MinValueValidator(0)], null=True, verbose_name='pK<sub>a</sub>', blank=True)),
                ('pi', models.FloatField(validators=[django.core.validators.MinValueValidator(0)], null=True, verbose_name='pI', blank=True)),
                ('log_p', models.FloatField(null=True, verbose_name='logP', blank=True)),
                ('log_d', models.FloatField(null=True, verbose_name='logD (pH 7.5)', blank=True)),
            ],
            options={
                'verbose_name_plural': 'Metabolites',
                'field_list': ['id', 'wid', 'name', 'traditional_name', 'iupac_name', 'synonyms', 'cross_references', 'type', 'empirical_formula', 'smiles', 'charge', 'is_hydrophobic', 'volume', 'deltag_formation', 'pka', 'pi', 'log_p', 'log_d', 'biomass_composition', 'media_composition', 'map_coordinates', 'comments', 'publication_references'],
                'verbose_name': 'Metabolite',
                'concrete_entry_model': True,
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'traditional_name', 'iupac_name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Structure', {'fields': [{'verbose_name': 'Structure', 'name': 'structure'}, 'empirical_formula', 'smiles', 'charge', 'is_hydrophobic', {'verbose_name': 'Molecular weight (Da)', 'name': 'molecular_weight'}, 'volume', 'deltag_formation', 'pka', 'pi', 'log_p', 'log_d']}), ('Concentrations', {'fields': ['biomass_composition', 'media_composition']}), ('Function', {'fields': [{'verbose_name': 'Coenzyme', 'name': 'coenzyme_participants'}, {'verbose_name': 'Complex subunit', 'name': 'protein_complex_biosythesis_participants'}, {'verbose_name': 'Prosthetic group', 'name': 'prosthetic_group_participants'}, {'verbose_name': 'Reaction participant', 'name': 'reaction_stoichiometry_participants'}]}), ('Parameters', {'fields': ['parameters']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
                'facet_fields': ['type', 'charge', 'is_hydrophobic'],
                'wid_unique': False,
            },
            bases=('cyano.molecule',),
        ),
        migrations.CreateModel(
            name='Protein',
            fields=[
                ('parent_ptr_molecule', models.OneToOneField(to='cyano.Molecule', verbose_name='Molecule', serialize=False, related_name='child_ptr_protein', parent_link=True, primary_key=True)),
            ],
            options={
                'verbose_name_plural': 'Proteins',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'prosthetic_groups', 'chaperones', 'dna_footprint', 'regulatory_rule', 'comments', 'publication_references'],
                'verbose_name': 'Protein',
                'concrete_entry_model': False,
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Structure', {'fields': ['prosthetic_groups', 'chaperones', 'dna_footprint']}), ('Regulation', {'fields': ['regulatory_rule']}), ('Function', {'fields': [{'verbose_name': 'Enzyme', 'name': 'enzyme_participants'}, {'verbose_name': 'Transcriptional regulation', 'name': 'transcriptional_regulations'}, {'verbose_name': 'Protein folding substrates', 'name': 'chaperone_substrates'}, {'verbose_name': 'Reaction participant', 'name': 'reaction_stoichiometry_participants'}, {'verbose_name': 'Complex subunit', 'name': 'protein_complex_biosythesis_participants'}]}), ('Parameters', {'fields': ['parameters']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
                'facet_fields': ['type', 'chaperones', 'dna_footprint__binding', 'dna_footprint__region'],
                'wid_unique': False,
            },
            bases=('cyano.molecule',),
        ),
        migrations.CreateModel(
            name='Stimulus',
            fields=[
                ('parent_ptr_molecule', models.OneToOneField(to='cyano.Molecule', verbose_name='Molecule', serialize=False, related_name='child_ptr_stimulus', parent_link=True, primary_key=True)),
                ('value', cyano.history.HistoryForeignKey(verbose_name='Value', related_name='+', to='cyano.EntryFloatData')),
            ],
            options={
                'verbose_name_plural': 'Stimuli',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'value', 'comments', 'publication_references'],
                'verbose_name': 'Stimulus',
                'concrete_entry_model': True,
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Value', {'fields': ['value']}), ('Function', {'fields': [{'verbose_name': 'Reaction participant', 'name': 'reaction_stoichiometry_participants'}, {'verbose_name': 'Complex subunit', 'name': 'protein_complex_biosythesis_participants'}]}), ('Parameters', {'fields': ['parameters']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
                'facet_fields': ['type', 'value__units'],
                'wid_unique': False,
            },
            bases=('cyano.molecule',),
        ),
        migrations.CreateModel(
            name='TranscriptionUnit',
            fields=[
                ('parent_ptr_molecule', models.OneToOneField(to='cyano.Molecule', verbose_name='Molecule', serialize=False, related_name='child_ptr_transcription_unit', parent_link=True, primary_key=True)),
                ('promoter_35_coordinate', models.IntegerField(null=True, verbose_name='Promoter -35 box coordinate (nt)', blank=True)),
                ('promoter_35_length', models.IntegerField(null=True, verbose_name='Promoter -35 box length (nt)', blank=True)),
                ('promoter_10_coordinate', models.IntegerField(null=True, verbose_name='Promoter -10 box coordinate (nt)', blank=True)),
                ('promoter_10_length', models.IntegerField(null=True, verbose_name='Promoter -10 box length (nt)', blank=True)),
                ('tss_coordinate', models.IntegerField(null=True, verbose_name='Transcription start site coordinate (nt)', blank=True)),
                ('genes', cyano.history.HistoryManyToManyField(verbose_name='Genes', related_name='transcription_units', to='cyano.Gene')),
            ],
            options={
                'verbose_name_plural': 'Transcription units',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'genes', 'promoter_35_coordinate', 'promoter_35_length', 'promoter_10_coordinate', 'promoter_10_length', 'tss_coordinate', 'comments', 'publication_references'],
                'verbose_name': 'Transcription unit',
                'concrete_entry_model': True,
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Structure (Hayflick media, 37C)', {'fields': [{'verbose_name': 'Structure', 'name': 'structure'}, {'verbose_name': 'Structure Filter', 'name': 'structure_filter'}, 'genes', 'promoter_35_coordinate', 'promoter_35_length', 'promoter_10_coordinate', 'promoter_10_length', 'tss_coordinate', {'verbose_name': 'Sequence', 'name': 'sequence'}]}), ('Regulation', {'fields': [{'verbose_name': 'Regulation', 'name': 'transcriptional_regulations'}]}), ('Function', {'fields': [{'verbose_name': 'Reaction participant', 'name': 'reaction_stoichiometry_participants'}, {'verbose_name': 'Complex subunit', 'name': 'protein_complex_biosythesis_participants'}]}), ('Parameters', {'fields': ['parameters']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
                'facet_fields': ['type'],
                'wid_unique': False,
            },
            bases=('cyano.molecule',),
        ),
        migrations.AddField(
            model_name='speciescomponent',
            name='publication_references',
            field=cyano.history.HistoryManyToManyField(verbose_name='Publications', blank=True, related_name='publication_referenced_components', to='cyano.PublicationReference'),
        ),
        migrations.AddField(
            model_name='speciescomponent',
            name='type',
            field=cyano.history.HistoryManyToManyField(verbose_name='Type', blank=True, related_name='members', to='cyano.Type'),
        ),
        migrations.AddField(
            model_name='species',
            name='publication_references',
            field=cyano.history.HistoryManyToManyField(verbose_name='Publication references', blank=True, related_name='publication_referenced_species', to='cyano.PublicationReference'),
        ),
        migrations.AddField(
            model_name='reactionstoichiometryparticipant',
            name='compartment',
            field=cyano.history.HistoryForeignKey(verbose_name='Compartment', related_name='+', to='cyano.Compartment'),
        ),
        migrations.AddField(
            model_name='reactionstoichiometryparticipant',
            name='molecule',
            field=cyano.history.HistoryForeignKey(verbose_name='Molecule', related_name='reaction_stoichiometry_participants', to='cyano.Molecule'),
        ),
        migrations.AddField(
            model_name='reaction',
            name='states',
            field=cyano.history.HistoryForeignKey(null=True, verbose_name='State', blank=True, related_name='reactions', to='cyano.State', on_delete=django.db.models.deletion.SET_NULL),
        ),
        migrations.AddField(
            model_name='reaction',
            name='stoichiometry',
            field=cyano.history.HistoryManyToManyField(verbose_name='Stoichiometry', related_name='reactions', to='cyano.ReactionStoichiometryParticipant'),
        ),
        migrations.AddField(
            model_name='proteincomplexbiosythesisparticipant',
            name='compartment',
            field=cyano.history.HistoryForeignKey(verbose_name='Compartment', related_name='+', to='cyano.Compartment'),
        ),
        migrations.AddField(
            model_name='proteincomplexbiosythesisparticipant',
            name='molecule',
            field=cyano.history.HistoryForeignKey(verbose_name='Molecule', related_name='protein_complex_biosythesis_participants', to='cyano.Molecule'),
        ),
        migrations.AddField(
            model_name='prostheticgroupparticipant',
            name='compartment',
            field=cyano.history.HistoryForeignKey(verbose_name='Compartment', related_name='+', to='cyano.Compartment'),
        ),
        migrations.AddField(
            model_name='parameter',
            name='molecules',
            field=cyano.history.HistoryManyToManyField(verbose_name='Molecules', blank=True, related_name='parameters', to='cyano.Molecule'),
        ),
        migrations.AddField(
            model_name='parameter',
            name='process',
            field=cyano.history.HistoryForeignKey(null=True, verbose_name='Process', blank=True, related_name='parameters', to='cyano.Process', on_delete=django.db.models.deletion.SET_NULL),
        ),
        migrations.AddField(
            model_name='parameter',
            name='reactions',
            field=cyano.history.HistoryManyToManyField(verbose_name='Reactions', blank=True, related_name='parameters', to='cyano.Reaction'),
        ),
        migrations.AddField(
            model_name='parameter',
            name='state',
            field=cyano.history.HistoryForeignKey(null=True, verbose_name='State', blank=True, related_name='parameters', to='cyano.State', on_delete=django.db.models.deletion.SET_NULL),
        ),
        migrations.AddField(
            model_name='parameter',
            name='value',
            field=cyano.history.HistoryForeignKey(to='cyano.EntryCharData', verbose_name='Value'),
        ),
        migrations.AddField(
            model_name='modificationreaction',
            name='compartment',
            field=cyano.history.HistoryForeignKey(verbose_name='Compartment', related_name='+', to='cyano.Compartment'),
        ),
        migrations.AddField(
            model_name='modificationreaction',
            name='molecule',
            field=cyano.history.HistoryForeignKey(verbose_name='Molecule', related_name='modification_reactions', to='cyano.Molecule'),
        ),
        migrations.AddField(
            model_name='metabolitemapcoordinate',
            name='compartment',
            field=cyano.history.HistoryForeignKey(verbose_name='Compartment', related_name='+', to='cyano.Compartment'),
        ),
        migrations.AddField(
            model_name='featureposition',
            name='chromosome_feature',
            field=cyano.history.HistoryForeignKey(verbose_name='', related_name='positions', to='cyano.ChromosomeFeature'),
        ),
        migrations.AddField(
            model_name='evidence',
            name='references',
            field=cyano.history.HistoryForeignKey(null=True, verbose_name='References', blank=True, related_name='evidence', to='cyano.PublicationReference'),
        ),
        migrations.AddField(
            model_name='enzymeparticipant',
            name='compartment',
            field=cyano.history.HistoryForeignKey(verbose_name='Compartment', related_name='+', to='cyano.Compartment'),
        ),
        migrations.AddField(
            model_name='coenzymeparticipant',
            name='compartment',
            field=cyano.history.HistoryForeignKey(verbose_name='Compartment', related_name='+', to='cyano.Compartment'),
        ),
        migrations.AddField(
            model_name='biomasscomposition',
            name='compartment',
            field=cyano.history.HistoryForeignKey(verbose_name='Compartment', related_name='biomass_compositions', to='cyano.Compartment'),
        ),
        migrations.CreateModel(
            name='Chromosome',
            fields=[
                ('parent_ptr_genome', models.OneToOneField(to='cyano.Genome', verbose_name='Genome', serialize=False, related_name='child_ptr_chromosome', parent_link=True, primary_key=True)),
            ],
            options={
                'verbose_name_plural': 'Chromosomes',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'sequence', 'length', 'comments', 'publication_references'],
                'verbose_name': 'Chromosome',
                'concrete_entry_model': True,
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Sequence', {'fields': [{'verbose_name': 'Structure Filter', 'name': 'structure_filter'}, {'verbose_name': 'Structure', 'name': 'structure'}, {'verbose_name': 'Extinction coefficient <br/>(260 nm, 25C, pH 7.0)', 'name': 'extinction_coefficient'}, {'verbose_name': 'pI', 'name': 'pi'}]}), ('Features', {'fields': ['genes', {'verbose_name': 'Other features', 'name': 'features'}]}), ('Function', {'fields': [{'verbose_name': 'Reaction participant', 'name': 'reaction_stoichiometry_participants'}, {'verbose_name': 'Complex subunit', 'name': 'protein_complex_biosythesis_participants'}]}), ('Parameters', {'fields': ['parameters']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
                'facet_fields': ['type'],
                'wid_unique': False,
            },
            bases=('cyano.genome',),
        ),
        migrations.CreateModel(
            name='MassSpectrometryProtein',
            fields=[
                ('parent_ptr_protein', models.OneToOneField(to='cyano.Protein', verbose_name='Protein', serialize=False, related_name='child_ptr_ms_protein', parent_link=True, primary_key=True)),
                ('score', models.FloatField(verbose_name='Protein Score')),
                ('coverage', models.FloatField(verbose_name='% Coverage')),
                ('sequence', models.TextField(verbose_name='Sequence', validators=[cyano.models.validate_protein_sequence])),
                ('length', models.PositiveIntegerField(verbose_name='Length (nt)')),
                ('pi', models.FloatField(verbose_name='Protein PI')),
                ('mass', models.FloatField(verbose_name='Protein Mass (Da)')),
                ('ambiguous', cyano.history.HistoryManyToManyField(verbose_name='Ambiguous Proteins', related_name='ambiguous', to='cyano.EntryBasicTextData')),
                ('sub', cyano.history.HistoryManyToManyField(verbose_name='Sub-Proteins', related_name='sub', to='cyano.EntryBasicTextData')),
            ],
            options={
                'verbose_name_plural': 'Mass Spectrometry Proteins',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'prosthetic_groups', 'chaperones', 'dna_footprint', 'regulatory_rule', 'score', 'coverage', 'pi', 'mass', 'sequence', 'comments', 'publication_references'],
                'verbose_name': 'Mass Spectrometry Protein',
                'concrete_entry_model': True,
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Related', {'fields': ['parent', {'verbose_name': 'Ambiguous Proteins', 'name': 'ambiguous_proteins'}, {'verbose_name': 'Sun-Proteins', 'name': 'sub_proteins'}]}), ('Structure', {'fields': ['prosthetic_groups', 'chaperones', 'dna_footprint', {'verbose_name': 'Sequence', 'name': 'sequence'}, {'verbose_name': 'Structure', 'name': 'structure'}]}), ('Regulation', {'fields': ['regulatory_rule']}), ('Function', {'fields': [{'verbose_name': 'Enzyme', 'name': 'enzyme_participants'}, {'verbose_name': 'Transcriptional regulation', 'name': 'transcriptional_regulations'}, {'verbose_name': 'Protein folding substrates', 'name': 'chaperone_substrates'}, {'verbose_name': 'Reaction participant', 'name': 'reaction_stoichiometry_participants'}, {'verbose_name': 'Complex subunit', 'name': 'protein_complex_biosythesis_participants'}]}), ('Statistics', {'fields': ['score', 'coverage', 'pi', 'mass']}), {'inline': 'protein_details'}, ('Parameters', {'fields': ['parameters']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
                'facet_fields': ['type', 'chaperones', 'dna_footprint__binding', 'dna_footprint__region'],
                'wid_unique': False,
            },
            bases=('cyano.protein',),
        ),
        migrations.CreateModel(
            name='Peptide',
            fields=[
                ('parent_ptr_protein', models.OneToOneField(to='cyano.Protein', verbose_name='Species component', serialize=False, related_name='child_ptr_peptide', parent_link=True, primary_key=True)),
                ('sequence', models.TextField(validators=[cyano.models.validate_protein_sequence], default='', verbose_name='Sequence', blank=True)),
                ('length', models.PositiveIntegerField(verbose_name='Length (nt)')),
                ('proteotypic', models.NullBooleanField(verbose_name='Proteotypic')),
                ('charge', models.IntegerField(verbose_name='Charge')),
                ('mass', models.FloatField(verbose_name='m/z')),
                ('zscore', models.FloatField(verbose_name='z-score')),
                ('retention_time', models.FloatField(verbose_name='Retention Time')),
                ('proteins', cyano.history.HistoryManyToManyField(verbose_name='Proteins belonging to the Peptide', related_name='peptides', to='cyano.EntryBasicTextData')),
            ],
            options={
                'verbose_name_plural': 'Peptides',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'prosthetic_groups', 'chaperones', 'dna_footprint', 'regulatory_rule', 'comments', 'publication_references'],
                'verbose_name': 'Peptide',
                'concrete_entry_model': True,
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Related', {'fields': ['parent', {'verbose_name': 'Matched Proteins', 'name': 'matched_proteins'}]}), ('Structure', {'fields': ['prosthetic_groups', 'chaperones', 'dna_footprint', {'verbose_name': 'Sequence', 'name': 'sequence'}]}), ('Regulation', {'fields': ['regulatory_rule']}), ('Function', {'fields': [{'verbose_name': 'Enzyme', 'name': 'enzyme_participants'}, {'verbose_name': 'Transcriptional regulation', 'name': 'transcriptional_regulations'}, {'verbose_name': 'Protein folding substrates', 'name': 'chaperone_substrates'}, {'verbose_name': 'Reaction participant', 'name': 'reaction_stoichiometry_participants'}, {'verbose_name': 'Complex subunit', 'name': 'protein_complex_biosythesis_participants'}]}), ('Statistics', {'fields': ['proteotypic', 'charge', 'mass', 'zscore', 'retention_time']}), ('Parameters', {'fields': ['parameters']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
                'facet_fields': ['type', 'chaperones', 'dna_footprint__binding', 'dna_footprint__region'],
                'wid_unique': False,
            },
            bases=('cyano.protein',),
        ),
        migrations.CreateModel(
            name='Plasmid',
            fields=[
                ('parent_ptr_genome', models.OneToOneField(to='cyano.Genome', verbose_name='Genome', serialize=False, related_name='child_ptr_plasmid', parent_link=True, primary_key=True)),
            ],
            options={
                'verbose_name_plural': 'Plasmids',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'sequence', 'length', 'comments', 'publication_references'],
                'verbose_name': 'Plasmid',
                'concrete_entry_model': True,
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Sequence', {'fields': [{'verbose_name': 'Structure Filter', 'name': 'structure_filter'}, {'verbose_name': 'Structure', 'name': 'structure'}, {'verbose_name': 'Extinction coefficient <br/>(260 nm, 25C, pH 7.0)', 'name': 'extinction_coefficient'}, {'verbose_name': 'pI', 'name': 'pi'}]}), ('Features', {'fields': ['genes', {'verbose_name': 'Other features', 'name': 'features'}]}), ('Function', {'fields': [{'verbose_name': 'Reaction participant', 'name': 'reaction_stoichiometry_participants'}, {'verbose_name': 'Complex subunit', 'name': 'protein_complex_biosythesis_participants'}]}), ('Parameters', {'fields': ['parameters']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
                'facet_fields': ['type'],
                'wid_unique': False,
            },
            bases=('cyano.genome',),
        ),
        migrations.CreateModel(
            name='ProteinComplex',
            fields=[
                ('parent_ptr_protein', models.OneToOneField(to='cyano.Protein', verbose_name='Protein', serialize=False, related_name='child_ptr_protein_complex', parent_link=True, primary_key=True)),
                ('biosynthesis', cyano.history.HistoryManyToManyField(verbose_name='Biosynthesis', related_name='protein_complexes', to='cyano.ProteinComplexBiosythesisParticipant')),
                ('disulfide_bonds', cyano.history.HistoryManyToManyField(verbose_name='Disulfide bonds (pH 7.5)', blank=True, related_name='protein_complexes', to='cyano.DisulfideBond')),
                ('formation_process', cyano.history.HistoryForeignKey(null=True, verbose_name='Formation process', blank=True, related_name='formed_complexes', to='cyano.Process', on_delete=django.db.models.deletion.SET_NULL)),
            ],
            options={
                'verbose_name_plural': 'Protein complexes',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'biosynthesis', 'disulfide_bonds', 'prosthetic_groups', 'dna_footprint', 'formation_process', 'chaperones', 'regulatory_rule', 'comments', 'publication_references'],
                'verbose_name': 'Protein complex',
                'concrete_entry_model': True,
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Structure', {'fields': ['biosynthesis', {'verbose_name': 'No. subunits', 'name': 'num_subunits'}, 'disulfide_bonds', 'prosthetic_groups', 'dna_footprint', {'verbose_name': 'Empirical formula (pH 7.5)', 'name': 'empirical_formula'}, {'verbose_name': 'Molecular weight (pH 7.5; Da)', 'name': 'molecular_weight'}, {'verbose_name': 'Extinction coefficient <br/>(260 nm, 25C, pH 7.0)', 'name': 'extinction_coefficient'}, {'verbose_name': 'Half life (OD (600 nm) = 0.3, <br/>M9 media, 36C; min)', 'name': 'half_life'}]}), ('Synthesis', {'fields': ['formation_process', 'chaperones', {'verbose_name': 'Localization', 'name': 'localization'}]}), ('Regulation', {'fields': ['regulatory_rule']}), ('Function', {'fields': [{'verbose_name': 'Enzyme', 'name': 'enzyme_participants'}, {'verbose_name': 'Transcriptional regulation', 'name': 'transcriptional_regulations'}, {'verbose_name': 'Protein folding substrates', 'name': 'chaperone_substrates'}, {'verbose_name': 'Reaction participant', 'name': 'reaction_stoichiometry_participants'}, {'verbose_name': 'Complex subunit', 'name': 'protein_complex_biosythesis_participants'}]}), ('Parameters', {'fields': ['parameters']}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
                'facet_fields': ['type', 'dna_footprint__binding', 'dna_footprint__region', 'formation_process', 'chaperones'],
                'wid_unique': False,
            },
            bases=('cyano.protein',),
        ),
        migrations.CreateModel(
            name='ProteinMonomer',
            fields=[
                ('parent_ptr_protein', models.OneToOneField(to='cyano.Protein', verbose_name='Protein', serialize=False, related_name='child_ptr_protein_monomer', parent_link=True, primary_key=True)),
            ],
            options={
                'verbose_name_plural': 'Protein monomers',
                'field_list': ['id', 'wid', 'name', 'synonyms', 'cross_references', 'type', 'gene', 'is_n_terminal_methionine_cleaved', 'signal_sequence', 'prosthetic_groups', 'dna_footprint', 'localization', 'chaperones', 'regulatory_rule', 'comments', 'publication_references'],
                'verbose_name': 'Protein monomer',
                'concrete_entry_model': True,
                'fieldsets': [('Type', {'fields': ['model_type']}), ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), ('Classification', {'fields': ['type']}), ('Genetics', {'fields': ['gene']}), ('Structure', {'fields': [{'verbose_name': 'Sequence', 'name': 'sequence'}, 'is_n_terminal_methionine_cleaved', 'signal_sequence', 'prosthetic_groups', 'disulfide_bonds', 'dna_footprint', {'verbose_name': 'Empirical formula (pH 7.5)', 'name': 'empirical_formula'}, {'verbose_name': 'Molecular weight (pH 7.5; Da)', 'name': 'molecular_weight'}, {'verbose_name': 'Extinction coefficient <br/>(260 nm, 25C, pH 7.0)', 'name': 'extinction_coefficient'}, {'verbose_name': 'pI', 'name': 'pi'}, {'verbose_name': 'Instability index', 'name': 'instability'}, {'verbose_name': 'Is stable', 'name': 'is_stable'}, {'verbose_name': 'Aliphatic index', 'name': 'aliphatic'}, {'verbose_name': 'GRAVY (25C, pH 7.0)', 'name': 'gravy'}, {'verbose_name': 'Half life (OD (600 nm) = 0.3, <br/>M9 media, 36C; min)', 'name': 'half_life'}]}), ('Synthesis', {'fields': ['localization', 'chaperones']}), ('Regulation', {'fields': ['regulatory_rule']}), ('Function', {'fields': [{'verbose_name': 'Enzyme', 'name': 'enzyme_participants'}, {'verbose_name': 'Transcriptional regulation', 'name': 'transcriptional_regulations'}, {'verbose_name': 'Protein folding substrates', 'name': 'chaperone_substrates'}, {'verbose_name': 'Reaction participant', 'name': 'reaction_stoichiometry_participants'}, {'verbose_name': 'Complex subunit', 'name': 'protein_complex_biosythesis_participants'}]}), ('Parameters', {'fields': ['parameters']}), ('Interactions', {'fields': [{'verbose_name': 'Protein/Metabolite interactions', 'name': 'interactions'}]}), ('Comments', {'fields': ['comments', 'publication_references']}), ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]})],
                'facet_fields': ['type', 'is_n_terminal_methionine_cleaved__value', 'signal_sequence__type', 'signal_sequence__location', 'dna_footprint__binding', 'dna_footprint__region', 'localization', 'chaperones'],
                'wid_unique': False,
            },
            bases=('cyano.protein',),
        ),
        migrations.AddField(
            model_name='transcriptionalregulation',
            name='transcription_factor',
            field=cyano.history.HistoryForeignKey(verbose_name='Transcripton factor', related_name='transcriptional_regulations', to='cyano.Protein'),
        ),
        migrations.AddField(
            model_name='transcriptionalregulation',
            name='transcription_unit',
            field=cyano.history.HistoryForeignKey(verbose_name='Transcription unit', related_name='transcriptional_regulations', to='cyano.TranscriptionUnit'),
        ),
        migrations.AddField(
            model_name='protein',
            name='chaperones',
            field=cyano.history.HistoryManyToManyField(verbose_name='Chaperones', blank=True, related_name='chaperone_substrates', to='cyano.Protein'),
        ),
        migrations.AddField(
            model_name='protein',
            name='dna_footprint',
            field=cyano.history.HistoryForeignKey(null=True, verbose_name='DNA footprint', blank=True, related_name='proteins', to='cyano.DNAFootprint'),
        ),
        migrations.AddField(
            model_name='protein',
            name='prosthetic_groups',
            field=cyano.history.HistoryManyToManyField(verbose_name='Prosthetic groups', blank=True, related_name='proteins', to='cyano.ProstheticGroupParticipant'),
        ),
        migrations.AddField(
            model_name='protein',
            name='regulatory_rule',
            field=cyano.history.HistoryForeignKey(null=True, verbose_name='Regulatory rule', blank=True, related_name='+', to='cyano.EntryCharData', on_delete=django.db.models.deletion.SET_NULL),
        ),
        migrations.AddField(
            model_name='prostheticgroupparticipant',
            name='metabolite',
            field=cyano.history.HistoryForeignKey(verbose_name='Metabolite', related_name='prosthetic_group_participants', to='cyano.Metabolite'),
        ),
        migrations.AddField(
            model_name='metabolite',
            name='biomass_composition',
            field=cyano.history.HistoryManyToManyField(verbose_name='Biomass composition (SP4 media, <br/>5% CO<sub>2</sub>, 37C; mmol gDCW<sup>-1</sup>)', blank=True, related_name='metabolites', to='cyano.BiomassComposition'),
        ),
        migrations.AddField(
            model_name='metabolite',
            name='map_coordinates',
            field=cyano.history.HistoryManyToManyField(verbose_name='Map coordinates', blank=True, related_name='metabolites', to='cyano.MetaboliteMapCoordinate'),
        ),
        migrations.AddField(
            model_name='metabolite',
            name='media_composition',
            field=cyano.history.HistoryForeignKey(null=True, verbose_name='Media composition (SP4; mM)', blank=True, related_name='metabolites', to='cyano.MediaComposition', on_delete=django.db.models.deletion.SET_NULL),
        ),
        migrations.AddField(
            model_name='gene',
            name='amino_acid',
            field=cyano.history.HistoryForeignKey(null=True, verbose_name='Amino acid', blank=True, related_name='genes', to='cyano.Metabolite', on_delete=django.db.models.deletion.SET_NULL),
        ),
        migrations.AddField(
            model_name='gene',
            name='chromosome',
            field=cyano.history.HistoryForeignKey(verbose_name='Chromosome or Plasmid', related_name='genes', to='cyano.Genome'),
        ),
        migrations.AddField(
            model_name='gene',
            name='codons',
            field=cyano.history.HistoryManyToManyField(verbose_name='Codons', blank=True, related_name='genes', to='cyano.Codon'),
        ),
        migrations.AddField(
            model_name='gene',
            name='expression',
            field=cyano.history.HistoryForeignKey(null=True, verbose_name='Relative expression', blank=True, related_name='+', to='cyano.EntryPositiveFloatData'),
        ),
        migrations.AddField(
            model_name='gene',
            name='half_life',
            field=cyano.history.HistoryForeignKey(null=True, verbose_name='Half life', blank=True, related_name='+', to='cyano.EntryPositiveFloatData'),
        ),
        migrations.AddField(
            model_name='gene',
            name='homologs',
            field=cyano.history.HistoryManyToManyField(verbose_name='Homologs', blank=True, related_name='genes', to='cyano.Homolog'),
        ),
        migrations.AddField(
            model_name='gene',
            name='is_essential',
            field=cyano.history.HistoryForeignKey(null=True, verbose_name='Is essential', blank=True, related_name='+', to='cyano.EntryBooleanData'),
        ),
        migrations.AddField(
            model_name='featureposition',
            name='chromosome',
            field=cyano.history.HistoryForeignKey(verbose_name='Chromosome or Plasmid', related_name='features', to='cyano.Genome'),
        ),
        migrations.AddField(
            model_name='enzymeparticipant',
            name='protein',
            field=cyano.history.HistoryForeignKey(verbose_name='Protein', related_name='enzyme_participants', to='cyano.Protein'),
        ),
        migrations.AddField(
            model_name='coenzymeparticipant',
            name='metabolite',
            field=cyano.history.HistoryForeignKey(verbose_name='Metabolite', related_name='coenzyme_participants', to='cyano.Metabolite'),
        ),
        migrations.AddField(
            model_name='proteinmonomer',
            name='gene',
            field=cyano.history.HistoryForeignKey(verbose_name='Gene', related_name='protein_monomers', to='cyano.Gene'),
        ),
        migrations.AddField(
            model_name='proteinmonomer',
            name='is_n_terminal_methionine_cleaved',
            field=cyano.history.HistoryForeignKey(null=True, verbose_name='Is N-terminal methionine cleaved', related_name='+', to='cyano.EntryBooleanData'),
        ),
        migrations.AddField(
            model_name='proteinmonomer',
            name='localization',
            field=cyano.history.HistoryForeignKey(null=True, verbose_name='Localization', related_name='protein_monomers', to='cyano.Compartment'),
        ),
        migrations.AddField(
            model_name='proteinmonomer',
            name='signal_sequence',
            field=cyano.history.HistoryForeignKey(null=True, verbose_name='Sequence sequence', blank=True, related_name='protein_monomers', to='cyano.SignalSequence', on_delete=django.db.models.deletion.SET_NULL),
        ),
        migrations.AddField(
            model_name='massspectrometryproteindetail',
            name='protein',
            field=cyano.history.HistoryForeignKey(related_name='protein_details', to='cyano.MassSpectrometryProtein'),
        ),
        migrations.AddField(
            model_name='disulfidebond',
            name='protein_monomer',
            field=cyano.history.HistoryForeignKey(verbose_name='Protein monomer', related_name='disulfide_bonds', to='cyano.ProteinMonomer'),
        ),
    ]
