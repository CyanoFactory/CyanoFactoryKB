# -*- coding: utf-8 -*-
import datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'Permission'
        db.create_table(u'cyano_permission', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(default=u'', max_length=255)),
            ('description', self.gf('django.db.models.fields.CharField')(default=u'', max_length=255)),
        ))
        db.send_create_signal(u'cyano', ['Permission'])

        # Adding model 'UserPermission'
        db.create_table(u'cyano_userpermission', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('entry', self.gf('cyano.history.HistoryForeignKey')(related_name=u'user_permissions', to=orm['cyano.Entry'])),
            ('user', self.gf('cyano.history.HistoryForeignKey')(related_name=u'permissions', to=orm['cyano.UserProfile'])),
        ))
        db.send_create_signal(u'cyano', ['UserPermission'])

        # Adding unique constraint on 'UserPermission', fields ['entry', 'user']
        db.create_unique(u'cyano_userpermission', ['entry_id', 'user_id'])

        # Adding M2M table for field allow on 'UserPermission'
        m2m_table_name = db.shorten_name(u'cyano_userpermission_allow')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('userpermission', models.ForeignKey(orm[u'cyano.userpermission'], null=False)),
            ('permission', models.ForeignKey(orm[u'cyano.permission'], null=False))
        ))
        db.create_unique(m2m_table_name, ['userpermission_id', 'permission_id'])

        # Adding M2M table for field deny on 'UserPermission'
        m2m_table_name = db.shorten_name(u'cyano_userpermission_deny')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('userpermission', models.ForeignKey(orm[u'cyano.userpermission'], null=False)),
            ('permission', models.ForeignKey(orm[u'cyano.permission'], null=False))
        ))
        db.create_unique(m2m_table_name, ['userpermission_id', 'permission_id'])

        # Adding model 'GroupPermission'
        db.create_table(u'cyano_grouppermission', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('entry', self.gf('cyano.history.HistoryForeignKey')(related_name=u'group_permissions', to=orm['cyano.Entry'])),
            ('group', self.gf('cyano.history.HistoryForeignKey')(related_name=u'permissions', to=orm['cyano.GroupProfile'])),
        ))
        db.send_create_signal(u'cyano', ['GroupPermission'])

        # Adding unique constraint on 'GroupPermission', fields ['entry', 'group']
        db.create_unique(u'cyano_grouppermission', ['entry_id', 'group_id'])

        # Adding M2M table for field allow on 'GroupPermission'
        m2m_table_name = db.shorten_name(u'cyano_grouppermission_allow')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('grouppermission', models.ForeignKey(orm[u'cyano.grouppermission'], null=False)),
            ('permission', models.ForeignKey(orm[u'cyano.permission'], null=False))
        ))
        db.create_unique(m2m_table_name, ['grouppermission_id', 'permission_id'])

        # Adding M2M table for field deny on 'GroupPermission'
        m2m_table_name = db.shorten_name(u'cyano_grouppermission_deny')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('grouppermission', models.ForeignKey(orm[u'cyano.grouppermission'], null=False)),
            ('permission', models.ForeignKey(orm[u'cyano.permission'], null=False))
        ))
        db.create_unique(m2m_table_name, ['grouppermission_id', 'permission_id'])

        # Adding model 'GroupProfile'
        db.create_table(u'cyano_groupprofile', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('group', self.gf('django.db.models.fields.related.OneToOneField')(related_name=u'profile', unique=True, to=orm['auth.Group'])),
            ('description', self.gf('django.db.models.fields.CharField')(default=u'', max_length=255, blank=True)),
        ))
        db.send_create_signal(u'cyano', ['GroupProfile'])

        # Adding model 'UserProfile'
        db.create_table(u'cyano_userprofile', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('user', self.gf('django.db.models.fields.related.OneToOneField')(related_name=u'profile', unique=True, to=orm['auth.User'])),
            ('affiliation', self.gf('django.db.models.fields.CharField')(default=u'', max_length=255, blank=True)),
            ('website', self.gf('django.db.models.fields.URLField')(default=u'', max_length=255, blank=True)),
            ('phone', self.gf('django.db.models.fields.CharField')(default=u'', max_length=255, blank=True)),
            ('address', self.gf('django.db.models.fields.CharField')(default=u'', max_length=255, blank=True)),
            ('city', self.gf('django.db.models.fields.CharField')(default=u'', max_length=255, blank=True)),
            ('state', self.gf('django.db.models.fields.CharField')(default=u'', max_length=255, blank=True)),
            ('zip', self.gf('django.db.models.fields.CharField')(default=u'', max_length=255, blank=True)),
            ('country', self.gf('django.db.models.fields.CharField')(default=u'', max_length=255, blank=True)),
            ('force_password_change', self.gf('django.db.models.fields.BooleanField')(default=False)),
        ))
        db.send_create_signal(u'cyano', ['UserProfile'])

        # Adding model 'TableMeta'
        db.create_table(u'cyano_tablemeta', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('table_name', self.gf('django.db.models.fields.CharField')(unique=True, max_length=255)),
            ('model_name', self.gf('django.db.models.fields.CharField')(unique=True, max_length=255)),
        ))
        db.send_create_signal(u'cyano', ['TableMeta'])

        # Adding model 'TableMetaColumn'
        db.create_table(u'cyano_tablemetacolumn', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('table', self.gf('cyano.history.HistoryForeignKey')(related_name=u'columns', to=orm['cyano.TableMeta'])),
            ('column_name', self.gf('django.db.models.fields.CharField')(max_length=255)),
            ('column_id', self.gf('django.db.models.fields.IntegerField')()),
        ))
        db.send_create_signal(u'cyano', ['TableMetaColumn'])

        # Adding unique constraint on 'TableMetaColumn', fields ['table', 'column_name']
        db.create_unique(u'cyano_tablemetacolumn', ['table_id', 'column_name'])

        # Adding model 'TableMetaManyToMany'
        db.create_table(u'cyano_tablemetamanytomany', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('m2m_table', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', unique=True, to=orm['cyano.TableMeta'])),
            ('source_table', self.gf('cyano.history.HistoryForeignKey')(related_name=u'm2ms_source', to=orm['cyano.TableMeta'])),
            ('target_table', self.gf('cyano.history.HistoryForeignKey')(related_name=u'm2ms_target', to=orm['cyano.TableMeta'])),
        ))
        db.send_create_signal(u'cyano', ['TableMetaManyToMany'])

        # Adding unique constraint on 'TableMetaManyToMany', fields ['m2m_table', 'source_table', 'target_table']
        db.create_unique(u'cyano_tablemetamanytomany', ['m2m_table_id', 'source_table_id', 'target_table_id'])

        # Adding model 'RevisionDetail'
        db.create_table(u'cyano_revisiondetail', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('user', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', to=orm['cyano.UserProfile'])),
            ('date', self.gf('django.db.models.fields.DateTimeField')(default=datetime.datetime.now)),
            ('reason', self.gf('django.db.models.fields.CharField')(default=u'', max_length=255, blank=True)),
        ))
        db.send_create_signal(u'cyano', ['RevisionDetail'])

        # Adding model 'Revision'
        db.create_table(u'cyano_revision', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('current', self.gf('cyano.history.HistoryForeignKey')(related_name=u'revisions', to=orm['cyano.Entry'])),
            ('detail', self.gf('cyano.history.HistoryForeignKey')(related_name=u'revisions', null=True, to=orm['cyano.RevisionDetail'])),
            ('action', self.gf('django.db.models.fields.CharField')(max_length=1)),
            ('column', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', to=orm['cyano.TableMetaColumn'])),
            ('new_value', self.gf('django.db.models.fields.TextField')(null=True, blank=True)),
        ))
        db.send_create_signal(u'cyano', ['Revision'])

        # Adding model 'RevisionManyToMany'
        db.create_table(u'cyano_revisionmanytomany', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('current', self.gf('cyano.history.HistoryForeignKey')(related_name=u'revisions_m2m', to=orm['cyano.Entry'])),
            ('detail', self.gf('cyano.history.HistoryForeignKey')(related_name=u'revisions_m2m', null=True, to=orm['cyano.RevisionDetail'])),
            ('action', self.gf('django.db.models.fields.CharField')(max_length=1)),
            ('table', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', to=orm['cyano.TableMetaManyToMany'])),
            ('new_value', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
        ))
        db.send_create_signal(u'cyano', ['RevisionManyToMany'])

        # Adding model 'Evidence'
        db.create_table(u'cyano_evidence', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('value', self.gf('django.db.models.fields.TextField')(default=u'', blank=True)),
            ('units', self.gf('django.db.models.fields.CharField')(max_length=255, null=True, blank=True)),
            ('is_experimentally_constrained', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('species', self.gf('django.db.models.fields.CharField')(max_length=255, null=True, blank=True)),
            ('media', self.gf('django.db.models.fields.CharField')(max_length=255, null=True, blank=True)),
            ('pH', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('temperature', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('comments', self.gf('django.db.models.fields.TextField')(default=u'', blank=True)),
            ('references', self.gf('cyano.history.HistoryForeignKey')(blank=True, related_name=u'evidence', null=True, to=orm['cyano.PublicationReference'])),
        ))
        db.send_create_signal(u'cyano', ['Evidence'])

        # Adding M2M table for field species_component on 'Evidence'
        m2m_table_name = db.shorten_name(u'cyano_evidence_species_component')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('evidence', models.ForeignKey(orm[u'cyano.evidence'], null=False)),
            ('speciescomponent', models.ForeignKey(orm[u'cyano.speciescomponent'], null=False))
        ))
        db.create_unique(m2m_table_name, ['evidence_id', 'speciescomponent_id'])

        # Adding model 'EntryBooleanData'
        db.create_table(u'cyano_entrybooleandata', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('detail', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', to=orm['cyano.RevisionDetail'])),
            ('value', self.gf('django.db.models.fields.BooleanField')(default=False)),
        ))
        db.send_create_signal(u'cyano', ['EntryBooleanData'])

        # Adding M2M table for field evidence on 'EntryBooleanData'
        m2m_table_name = db.shorten_name(u'cyano_entrybooleandata_evidence')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('entrybooleandata', models.ForeignKey(orm[u'cyano.entrybooleandata'], null=False)),
            ('evidence', models.ForeignKey(orm[u'cyano.evidence'], null=False))
        ))
        db.create_unique(m2m_table_name, ['entrybooleandata_id', 'evidence_id'])

        # Adding model 'EntryCharData'
        db.create_table(u'cyano_entrychardata', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('detail', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', to=orm['cyano.RevisionDetail'])),
            ('value', self.gf('django.db.models.fields.CharField')(default=u'', max_length=255, blank=True)),
            ('units', self.gf('django.db.models.fields.CharField')(default=u'', max_length=255, blank=True)),
        ))
        db.send_create_signal(u'cyano', ['EntryCharData'])

        # Adding M2M table for field evidence on 'EntryCharData'
        m2m_table_name = db.shorten_name(u'cyano_entrychardata_evidence')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('entrychardata', models.ForeignKey(orm[u'cyano.entrychardata'], null=False)),
            ('evidence', models.ForeignKey(orm[u'cyano.evidence'], null=False))
        ))
        db.create_unique(m2m_table_name, ['entrychardata_id', 'evidence_id'])

        # Adding model 'EntryFloatData'
        db.create_table(u'cyano_entryfloatdata', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('detail', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', to=orm['cyano.RevisionDetail'])),
            ('value', self.gf('django.db.models.fields.FloatField')()),
            ('units', self.gf('django.db.models.fields.CharField')(default=u'', max_length=255, blank=True)),
        ))
        db.send_create_signal(u'cyano', ['EntryFloatData'])

        # Adding M2M table for field evidence on 'EntryFloatData'
        m2m_table_name = db.shorten_name(u'cyano_entryfloatdata_evidence')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('entryfloatdata', models.ForeignKey(orm[u'cyano.entryfloatdata'], null=False)),
            ('evidence', models.ForeignKey(orm[u'cyano.evidence'], null=False))
        ))
        db.create_unique(m2m_table_name, ['entryfloatdata_id', 'evidence_id'])

        # Adding model 'EntryPositiveFloatData'
        db.create_table(u'cyano_entrypositivefloatdata', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('detail', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', to=orm['cyano.RevisionDetail'])),
            ('value', self.gf('django.db.models.fields.FloatField')()),
            ('units', self.gf('django.db.models.fields.CharField')(default=u'', max_length=255, blank=True)),
        ))
        db.send_create_signal(u'cyano', ['EntryPositiveFloatData'])

        # Adding M2M table for field evidence on 'EntryPositiveFloatData'
        m2m_table_name = db.shorten_name(u'cyano_entrypositivefloatdata_evidence')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('entrypositivefloatdata', models.ForeignKey(orm[u'cyano.entrypositivefloatdata'], null=False)),
            ('evidence', models.ForeignKey(orm[u'cyano.evidence'], null=False))
        ))
        db.create_unique(m2m_table_name, ['entrypositivefloatdata_id', 'evidence_id'])

        # Adding model 'EntryTextData'
        db.create_table(u'cyano_entrytextdata', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('detail', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', to=orm['cyano.RevisionDetail'])),
            ('value', self.gf('django.db.models.fields.TextField')(default=u'', blank=True)),
            ('units', self.gf('django.db.models.fields.CharField')(default=u'', max_length=255, blank=True)),
        ))
        db.send_create_signal(u'cyano', ['EntryTextData'])

        # Adding M2M table for field evidence on 'EntryTextData'
        m2m_table_name = db.shorten_name(u'cyano_entrytextdata_evidence')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('entrytextdata', models.ForeignKey(orm[u'cyano.entrytextdata'], null=False)),
            ('evidence', models.ForeignKey(orm[u'cyano.evidence'], null=False))
        ))
        db.create_unique(m2m_table_name, ['entrytextdata_id', 'evidence_id'])

        # Adding model 'BindingSite'
        db.create_table(u'cyano_bindingsite', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('detail', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', to=orm['cyano.RevisionDetail'])),
            ('coordinate', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('length', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('direction', self.gf('django.db.models.fields.CharField')(max_length=10)),
        ))
        db.send_create_signal(u'cyano', ['BindingSite'])

        # Adding M2M table for field evidence on 'BindingSite'
        m2m_table_name = db.shorten_name(u'cyano_bindingsite_evidence')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('bindingsite', models.ForeignKey(orm[u'cyano.bindingsite'], null=False)),
            ('evidence', models.ForeignKey(orm[u'cyano.evidence'], null=False))
        ))
        db.create_unique(m2m_table_name, ['bindingsite_id', 'evidence_id'])

        # Adding model 'BiomassComposition'
        db.create_table(u'cyano_biomasscomposition', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('detail', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', to=orm['cyano.RevisionDetail'])),
            ('concentration', self.gf('django.db.models.fields.FloatField')()),
            ('compartment', self.gf('cyano.history.HistoryForeignKey')(related_name=u'biomass_compositions', to=orm['cyano.Compartment'])),
        ))
        db.send_create_signal(u'cyano', ['BiomassComposition'])

        # Adding M2M table for field evidence on 'BiomassComposition'
        m2m_table_name = db.shorten_name(u'cyano_biomasscomposition_evidence')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('biomasscomposition', models.ForeignKey(orm[u'cyano.biomasscomposition'], null=False)),
            ('evidence', models.ForeignKey(orm[u'cyano.evidence'], null=False))
        ))
        db.create_unique(m2m_table_name, ['biomasscomposition_id', 'evidence_id'])

        # Adding model 'Codon'
        db.create_table(u'cyano_codon', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('detail', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', to=orm['cyano.RevisionDetail'])),
            ('sequence', self.gf('django.db.models.fields.CharField')(max_length=3)),
        ))
        db.send_create_signal(u'cyano', ['Codon'])

        # Adding M2M table for field evidence on 'Codon'
        m2m_table_name = db.shorten_name(u'cyano_codon_evidence')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('codon', models.ForeignKey(orm[u'cyano.codon'], null=False)),
            ('evidence', models.ForeignKey(orm[u'cyano.evidence'], null=False))
        ))
        db.create_unique(m2m_table_name, ['codon_id', 'evidence_id'])

        # Adding model 'CoenzymeParticipant'
        db.create_table(u'cyano_coenzymeparticipant', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('detail', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', to=orm['cyano.RevisionDetail'])),
            ('metabolite', self.gf('cyano.history.HistoryForeignKey')(related_name=u'coenzyme_participants', to=orm['cyano.Metabolite'])),
            ('compartment', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', to=orm['cyano.Compartment'])),
            ('coefficient', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
        ))
        db.send_create_signal(u'cyano', ['CoenzymeParticipant'])

        # Adding M2M table for field evidence on 'CoenzymeParticipant'
        m2m_table_name = db.shorten_name(u'cyano_coenzymeparticipant_evidence')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('coenzymeparticipant', models.ForeignKey(orm[u'cyano.coenzymeparticipant'], null=False)),
            ('evidence', models.ForeignKey(orm[u'cyano.evidence'], null=False))
        ))
        db.create_unique(m2m_table_name, ['coenzymeparticipant_id', 'evidence_id'])

        # Adding model 'DisulfideBond'
        db.create_table(u'cyano_disulfidebond', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('detail', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', to=orm['cyano.RevisionDetail'])),
            ('protein_monomer', self.gf('cyano.history.HistoryForeignKey')(related_name=u'disulfide_bonds', to=orm['cyano.ProteinMonomer'])),
            ('residue_1', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('residue_2', self.gf('django.db.models.fields.PositiveIntegerField')()),
        ))
        db.send_create_signal(u'cyano', ['DisulfideBond'])

        # Adding M2M table for field evidence on 'DisulfideBond'
        m2m_table_name = db.shorten_name(u'cyano_disulfidebond_evidence')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('disulfidebond', models.ForeignKey(orm[u'cyano.disulfidebond'], null=False)),
            ('evidence', models.ForeignKey(orm[u'cyano.evidence'], null=False))
        ))
        db.create_unique(m2m_table_name, ['disulfidebond_id', 'evidence_id'])

        # Adding model 'DNAFootprint'
        db.create_table(u'cyano_dnafootprint', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('detail', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', to=orm['cyano.RevisionDetail'])),
            ('length', self.gf('django.db.models.fields.PositiveIntegerField')(null=True, blank=True)),
            ('binding', self.gf('django.db.models.fields.CharField')(default=u'', max_length=10, blank=True)),
            ('region', self.gf('django.db.models.fields.CharField')(default=u'', max_length=10, blank=True)),
        ))
        db.send_create_signal(u'cyano', ['DNAFootprint'])

        # Adding M2M table for field evidence on 'DNAFootprint'
        m2m_table_name = db.shorten_name(u'cyano_dnafootprint_evidence')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('dnafootprint', models.ForeignKey(orm[u'cyano.dnafootprint'], null=False)),
            ('evidence', models.ForeignKey(orm[u'cyano.evidence'], null=False))
        ))
        db.create_unique(m2m_table_name, ['dnafootprint_id', 'evidence_id'])

        # Adding model 'EnzymeParticipant'
        db.create_table(u'cyano_enzymeparticipant', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('detail', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', to=orm['cyano.RevisionDetail'])),
            ('protein', self.gf('cyano.history.HistoryForeignKey')(related_name=u'enzyme_participants', to=orm['cyano.Protein'])),
            ('compartment', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', to=orm['cyano.Compartment'])),
        ))
        db.send_create_signal(u'cyano', ['EnzymeParticipant'])

        # Adding M2M table for field evidence on 'EnzymeParticipant'
        m2m_table_name = db.shorten_name(u'cyano_enzymeparticipant_evidence')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('enzymeparticipant', models.ForeignKey(orm[u'cyano.enzymeparticipant'], null=False)),
            ('evidence', models.ForeignKey(orm[u'cyano.evidence'], null=False))
        ))
        db.create_unique(m2m_table_name, ['enzymeparticipant_id', 'evidence_id'])

        # Adding model 'Homolog'
        db.create_table(u'cyano_homolog', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('detail', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', to=orm['cyano.RevisionDetail'])),
            ('xid', self.gf('django.db.models.fields.CharField')(max_length=255)),
            ('species', self.gf('django.db.models.fields.CharField')(max_length=20)),
        ))
        db.send_create_signal(u'cyano', ['Homolog'])

        # Adding M2M table for field evidence on 'Homolog'
        m2m_table_name = db.shorten_name(u'cyano_homolog_evidence')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('homolog', models.ForeignKey(orm[u'cyano.homolog'], null=False)),
            ('evidence', models.ForeignKey(orm[u'cyano.evidence'], null=False))
        ))
        db.create_unique(m2m_table_name, ['homolog_id', 'evidence_id'])

        # Adding model 'Kinetics'
        db.create_table(u'cyano_kinetics', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('detail', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', to=orm['cyano.RevisionDetail'])),
            ('rate_law', self.gf('django.db.models.fields.CharField')(default=u'', max_length=255, blank=True)),
            ('km', self.gf('django.db.models.fields.CharField')(max_length=255, blank=True)),
            ('vmax', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('vmax_unit', self.gf('django.db.models.fields.CharField')(max_length=255, blank=True)),
        ))
        db.send_create_signal(u'cyano', ['Kinetics'])

        # Adding M2M table for field evidence on 'Kinetics'
        m2m_table_name = db.shorten_name(u'cyano_kinetics_evidence')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('kinetics', models.ForeignKey(orm[u'cyano.kinetics'], null=False)),
            ('evidence', models.ForeignKey(orm[u'cyano.evidence'], null=False))
        ))
        db.create_unique(m2m_table_name, ['kinetics_id', 'evidence_id'])

        # Adding model 'MediaComposition'
        db.create_table(u'cyano_mediacomposition', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('detail', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', to=orm['cyano.RevisionDetail'])),
            ('concentration', self.gf('django.db.models.fields.FloatField')()),
            ('is_diffused', self.gf('django.db.models.fields.BooleanField')(default=False)),
        ))
        db.send_create_signal(u'cyano', ['MediaComposition'])

        # Adding M2M table for field evidence on 'MediaComposition'
        m2m_table_name = db.shorten_name(u'cyano_mediacomposition_evidence')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('mediacomposition', models.ForeignKey(orm[u'cyano.mediacomposition'], null=False)),
            ('evidence', models.ForeignKey(orm[u'cyano.evidence'], null=False))
        ))
        db.create_unique(m2m_table_name, ['mediacomposition_id', 'evidence_id'])

        # Adding model 'MetaboliteMapCoordinate'
        db.create_table(u'cyano_metabolitemapcoordinate', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('compartment', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', to=orm['cyano.Compartment'])),
            ('x', self.gf('django.db.models.fields.FloatField')()),
            ('y', self.gf('django.db.models.fields.FloatField')()),
        ))
        db.send_create_signal(u'cyano', ['MetaboliteMapCoordinate'])

        # Adding model 'ModificationReaction'
        db.create_table(u'cyano_modificationreaction', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('detail', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', to=orm['cyano.RevisionDetail'])),
            ('molecule', self.gf('cyano.history.HistoryForeignKey')(related_name=u'modification_reactions', to=orm['cyano.Molecule'])),
            ('compartment', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', to=orm['cyano.Compartment'])),
            ('position', self.gf('django.db.models.fields.PositiveIntegerField')(null=True, blank=True)),
        ))
        db.send_create_signal(u'cyano', ['ModificationReaction'])

        # Adding M2M table for field evidence on 'ModificationReaction'
        m2m_table_name = db.shorten_name(u'cyano_modificationreaction_evidence')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('modificationreaction', models.ForeignKey(orm[u'cyano.modificationreaction'], null=False)),
            ('evidence', models.ForeignKey(orm[u'cyano.evidence'], null=False))
        ))
        db.create_unique(m2m_table_name, ['modificationreaction_id', 'evidence_id'])

        # Adding model 'ProstheticGroupParticipant'
        db.create_table(u'cyano_prostheticgroupparticipant', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('detail', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', to=orm['cyano.RevisionDetail'])),
            ('metabolite', self.gf('cyano.history.HistoryForeignKey')(related_name=u'prosthetic_group_participants', to=orm['cyano.Metabolite'])),
            ('compartment', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', to=orm['cyano.Compartment'])),
            ('coefficient', self.gf('django.db.models.fields.PositiveIntegerField')(null=True, blank=True)),
        ))
        db.send_create_signal(u'cyano', ['ProstheticGroupParticipant'])

        # Adding M2M table for field evidence on 'ProstheticGroupParticipant'
        m2m_table_name = db.shorten_name(u'cyano_prostheticgroupparticipant_evidence')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('prostheticgroupparticipant', models.ForeignKey(orm[u'cyano.prostheticgroupparticipant'], null=False)),
            ('evidence', models.ForeignKey(orm[u'cyano.evidence'], null=False))
        ))
        db.create_unique(m2m_table_name, ['prostheticgroupparticipant_id', 'evidence_id'])

        # Adding model 'ProteinComplexBiosythesisParticipant'
        db.create_table(u'cyano_proteincomplexbiosythesisparticipant', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('detail', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', to=orm['cyano.RevisionDetail'])),
            ('molecule', self.gf('cyano.history.HistoryForeignKey')(related_name=u'protein_complex_biosythesis_participants', to=orm['cyano.Molecule'])),
            ('residue', self.gf('django.db.models.fields.PositiveIntegerField')(null=True, blank=True)),
            ('coefficient', self.gf('django.db.models.fields.FloatField')()),
            ('compartment', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', to=orm['cyano.Compartment'])),
        ))
        db.send_create_signal(u'cyano', ['ProteinComplexBiosythesisParticipant'])

        # Adding M2M table for field evidence on 'ProteinComplexBiosythesisParticipant'
        m2m_table_name = db.shorten_name(u'cyano_proteincomplexbiosythesisparticipant_evidence')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('proteincomplexbiosythesisparticipant', models.ForeignKey(orm[u'cyano.proteincomplexbiosythesisparticipant'], null=False)),
            ('evidence', models.ForeignKey(orm[u'cyano.evidence'], null=False))
        ))
        db.create_unique(m2m_table_name, ['proteincomplexbiosythesisparticipant_id', 'evidence_id'])

        # Adding model 'ReactionMapCoordinate'
        db.create_table(u'cyano_reactionmapcoordinate', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('path', self.gf('django.db.models.fields.TextField')()),
            ('value_x', self.gf('django.db.models.fields.FloatField')()),
            ('value_y', self.gf('django.db.models.fields.FloatField')()),
            ('label_x', self.gf('django.db.models.fields.FloatField')()),
            ('label_y', self.gf('django.db.models.fields.FloatField')()),
        ))
        db.send_create_signal(u'cyano', ['ReactionMapCoordinate'])

        # Adding model 'ReactionStoichiometryParticipant'
        db.create_table(u'cyano_reactionstoichiometryparticipant', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('detail', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', to=orm['cyano.RevisionDetail'])),
            ('molecule', self.gf('cyano.history.HistoryForeignKey')(related_name=u'reaction_stoichiometry_participants', to=orm['cyano.Molecule'])),
            ('coefficient', self.gf('django.db.models.fields.FloatField')()),
            ('compartment', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', to=orm['cyano.Compartment'])),
        ))
        db.send_create_signal(u'cyano', ['ReactionStoichiometryParticipant'])

        # Adding M2M table for field evidence on 'ReactionStoichiometryParticipant'
        m2m_table_name = db.shorten_name(u'cyano_reactionstoichiometryparticipant_evidence')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('reactionstoichiometryparticipant', models.ForeignKey(orm[u'cyano.reactionstoichiometryparticipant'], null=False)),
            ('evidence', models.ForeignKey(orm[u'cyano.evidence'], null=False))
        ))
        db.create_unique(m2m_table_name, ['reactionstoichiometryparticipant_id', 'evidence_id'])

        # Adding model 'SignalSequence'
        db.create_table(u'cyano_signalsequence', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('detail', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', to=orm['cyano.RevisionDetail'])),
            ('type', self.gf('django.db.models.fields.CharField')(max_length=20)),
            ('location', self.gf('django.db.models.fields.CharField')(max_length=1)),
            ('length', self.gf('django.db.models.fields.PositiveIntegerField')()),
        ))
        db.send_create_signal(u'cyano', ['SignalSequence'])

        # Adding M2M table for field evidence on 'SignalSequence'
        m2m_table_name = db.shorten_name(u'cyano_signalsequence_evidence')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('signalsequence', models.ForeignKey(orm[u'cyano.signalsequence'], null=False)),
            ('evidence', models.ForeignKey(orm[u'cyano.evidence'], null=False))
        ))
        db.create_unique(m2m_table_name, ['signalsequence_id', 'evidence_id'])

        # Adding model 'Synonym'
        db.create_table(u'cyano_synonym', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=255)),
        ))
        db.send_create_signal(u'cyano', ['Synonym'])

        # Adding model 'Entry'
        db.create_table(u'cyano_entry', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('model_type', self.gf('cyano.history.HistoryForeignKey')(to=orm['cyano.TableMeta'])),
            ('wid', self.gf('django.db.models.fields.SlugField')(max_length=150)),
            ('name', self.gf('django.db.models.fields.CharField')(default=u'', max_length=255, blank=True)),
            ('comments', self.gf('django.db.models.fields.TextField')(default=u'', blank=True)),
            ('created_detail', self.gf('cyano.history.HistoryForeignKey')(related_name=u'entry_created_detail', to=orm['cyano.RevisionDetail'])),
            ('detail', self.gf('cyano.history.HistoryForeignKey')(related_name=u'entry_detail', to=orm['cyano.RevisionDetail'])),
        ))
        db.send_create_signal(u'cyano', ['Entry'])

        # Adding M2M table for field synonyms on 'Entry'
        m2m_table_name = db.shorten_name(u'cyano_entry_synonyms')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('entry', models.ForeignKey(orm[u'cyano.entry'], null=False)),
            ('synonym', models.ForeignKey(orm[u'cyano.synonym'], null=False))
        ))
        db.create_unique(m2m_table_name, ['entry_id', 'synonym_id'])

        # Adding model 'CrossReferenceMeta'
        db.create_table(u'cyano_crossreferencemeta', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=255)),
        ))
        db.send_create_signal(u'cyano', ['CrossReferenceMeta'])

        # Adding model 'SpeciesComponent'
        db.create_table(u'cyano_speciescomponent', (
            (u'entry_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['cyano.Entry'], unique=True, primary_key=True)),
        ))
        db.send_create_signal(u'cyano', ['SpeciesComponent'])

        # Adding M2M table for field species on 'SpeciesComponent'
        m2m_table_name = db.shorten_name(u'cyano_speciescomponent_species')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('speciescomponent', models.ForeignKey(orm[u'cyano.speciescomponent'], null=False)),
            ('species', models.ForeignKey(orm[u'cyano.species'], null=False))
        ))
        db.create_unique(m2m_table_name, ['speciescomponent_id', 'species_id'])

        # Adding M2M table for field type on 'SpeciesComponent'
        m2m_table_name = db.shorten_name(u'cyano_speciescomponent_type')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('speciescomponent', models.ForeignKey(orm[u'cyano.speciescomponent'], null=False)),
            ('type', models.ForeignKey(orm[u'cyano.type'], null=False))
        ))
        db.create_unique(m2m_table_name, ['speciescomponent_id', 'type_id'])

        # Adding M2M table for field cross_references on 'SpeciesComponent'
        m2m_table_name = db.shorten_name(u'cyano_speciescomponent_cross_references')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('speciescomponent', models.ForeignKey(orm[u'cyano.speciescomponent'], null=False)),
            ('crossreference', models.ForeignKey(orm[u'cyano.crossreference'], null=False))
        ))
        db.create_unique(m2m_table_name, ['speciescomponent_id', 'crossreference_id'])

        # Adding M2M table for field publication_references on 'SpeciesComponent'
        m2m_table_name = db.shorten_name(u'cyano_speciescomponent_publication_references')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('speciescomponent', models.ForeignKey(orm[u'cyano.speciescomponent'], null=False)),
            ('publicationreference', models.ForeignKey(orm[u'cyano.publicationreference'], null=False))
        ))
        db.create_unique(m2m_table_name, ['speciescomponent_id', 'publicationreference_id'])

        # Adding model 'Molecule'
        db.create_table(u'cyano_molecule', (
            ('parent_ptr_species_component', self.gf('django.db.models.fields.related.OneToOneField')(related_name=u'child_ptr_molecule', unique=True, primary_key=True, to=orm['cyano.SpeciesComponent'])),
        ))
        db.send_create_signal(u'cyano', ['Molecule'])

        # Adding model 'Protein'
        db.create_table(u'cyano_protein', (
            ('parent_ptr_molecule', self.gf('django.db.models.fields.related.OneToOneField')(related_name=u'child_ptr_protein', unique=True, primary_key=True, to=orm['cyano.Molecule'])),
            ('dna_footprint', self.gf('cyano.history.HistoryForeignKey')(blank=True, related_name=u'proteins', null=True, to=orm['cyano.DNAFootprint'])),
            ('regulatory_rule', self.gf('cyano.history.HistoryForeignKey')(blank=True, related_name=u'+', null=True, on_delete=models.SET_NULL, to=orm['cyano.EntryCharData'])),
        ))
        db.send_create_signal(u'cyano', ['Protein'])

        # Adding M2M table for field prosthetic_groups on 'Protein'
        m2m_table_name = db.shorten_name(u'cyano_protein_prosthetic_groups')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('protein', models.ForeignKey(orm[u'cyano.protein'], null=False)),
            ('prostheticgroupparticipant', models.ForeignKey(orm[u'cyano.prostheticgroupparticipant'], null=False))
        ))
        db.create_unique(m2m_table_name, ['protein_id', 'prostheticgroupparticipant_id'])

        # Adding M2M table for field chaperones on 'Protein'
        m2m_table_name = db.shorten_name(u'cyano_protein_chaperones')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('from_protein', models.ForeignKey(orm[u'cyano.protein'], null=False)),
            ('to_protein', models.ForeignKey(orm[u'cyano.protein'], null=False))
        ))
        db.create_unique(m2m_table_name, ['from_protein_id', 'to_protein_id'])

        # Adding model 'DNA'
        db.create_table(u'cyano_dna', (
            ('parent_ptr_molecule', self.gf('django.db.models.fields.related.OneToOneField')(related_name=u'child_ptr_dna', unique=True, primary_key=True, to=orm['cyano.Molecule'])),
            ('sequence', self.gf('django.db.models.fields.TextField')(default=u'', blank=True)),
            ('length', self.gf('django.db.models.fields.PositiveIntegerField')()),
        ))
        db.send_create_signal(u'cyano', ['DNA'])

        # Adding model 'Chromosome'
        db.create_table(u'cyano_chromosome', (
            ('parent_ptr_dna', self.gf('django.db.models.fields.related.OneToOneField')(related_name=u'child_ptr_chromosome', unique=True, primary_key=True, to=orm['cyano.DNA'])),
        ))
        db.send_create_signal(u'cyano', ['Chromosome'])

        # Adding model 'Plasmid'
        db.create_table(u'cyano_plasmid', (
            ('parent_ptr_dna', self.gf('django.db.models.fields.related.OneToOneField')(related_name=u'child_ptr_plasmid', unique=True, primary_key=True, to=orm['cyano.DNA'])),
        ))
        db.send_create_signal(u'cyano', ['Plasmid'])

        # Adding model 'ChromosomeFeature'
        db.create_table(u'cyano_chromosomefeature', (
            ('parent_ptr_species_component', self.gf('django.db.models.fields.related.OneToOneField')(related_name=u'child_ptr_chromosome_feature', unique=True, primary_key=True, to=orm['cyano.SpeciesComponent'])),
            ('chromosome', self.gf('cyano.history.HistoryForeignKey')(related_name=u'features', to=orm['cyano.DNA'])),
            ('coordinate', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('length', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('direction', self.gf('django.db.models.fields.CharField')(max_length=10)),
        ))
        db.send_create_signal(u'cyano', ['ChromosomeFeature'])

        # Adding model 'Compartment'
        db.create_table(u'cyano_compartment', (
            ('parent_ptr_species_component', self.gf('django.db.models.fields.related.OneToOneField')(related_name=u'child_ptr_compartment', unique=True, primary_key=True, to=orm['cyano.SpeciesComponent'])),
        ))
        db.send_create_signal(u'cyano', ['Compartment'])

        # Adding model 'Gene'
        db.create_table(u'cyano_gene', (
            ('parent_ptr_molecule', self.gf('django.db.models.fields.related.OneToOneField')(related_name=u'child_ptr_gene', unique=True, primary_key=True, to=orm['cyano.Molecule'])),
            ('symbol', self.gf('django.db.models.fields.CharField')(default=u'', max_length=255, blank=True)),
            ('chromosome', self.gf('cyano.history.HistoryForeignKey')(related_name=u'genes', to=orm['cyano.DNA'])),
            ('coordinate', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('length', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('direction', self.gf('django.db.models.fields.CharField')(max_length=10)),
            ('is_essential', self.gf('cyano.history.HistoryForeignKey')(blank=True, related_name=u'+', null=True, to=orm['cyano.EntryBooleanData'])),
            ('expression', self.gf('cyano.history.HistoryForeignKey')(blank=True, related_name=u'+', null=True, to=orm['cyano.EntryPositiveFloatData'])),
            ('half_life', self.gf('cyano.history.HistoryForeignKey')(blank=True, related_name=u'+', null=True, to=orm['cyano.EntryPositiveFloatData'])),
            ('amino_acid', self.gf('cyano.history.HistoryForeignKey')(blank=True, related_name=u'genes', null=True, on_delete=models.SET_NULL, to=orm['cyano.Metabolite'])),
        ))
        db.send_create_signal(u'cyano', ['Gene'])

        # Adding M2M table for field codons on 'Gene'
        m2m_table_name = db.shorten_name(u'cyano_gene_codons')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('gene', models.ForeignKey(orm[u'cyano.gene'], null=False)),
            ('codon', models.ForeignKey(orm[u'cyano.codon'], null=False))
        ))
        db.create_unique(m2m_table_name, ['gene_id', 'codon_id'])

        # Adding M2M table for field homologs on 'Gene'
        m2m_table_name = db.shorten_name(u'cyano_gene_homologs')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('gene', models.ForeignKey(orm[u'cyano.gene'], null=False)),
            ('homolog', models.ForeignKey(orm[u'cyano.homolog'], null=False))
        ))
        db.create_unique(m2m_table_name, ['gene_id', 'homolog_id'])

        # Adding model 'Metabolite'
        db.create_table(u'cyano_metabolite', (
            ('parent_ptr_molecule', self.gf('django.db.models.fields.related.OneToOneField')(related_name=u'child_ptr_metabolite', unique=True, primary_key=True, to=orm['cyano.Molecule'])),
            ('traditional_name', self.gf('django.db.models.fields.TextField')(default=u'', blank=True)),
            ('iupac_name', self.gf('django.db.models.fields.TextField')(default=u'', blank=True)),
            ('empirical_formula', self.gf('django.db.models.fields.TextField')()),
            ('smiles', self.gf('django.db.models.fields.TextField')(default=u'', blank=True)),
            ('charge', self.gf('django.db.models.fields.IntegerField')()),
            ('is_hydrophobic', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('volume', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('deltag_formation', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('pka', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('pi', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('log_p', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('log_d', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('media_composition', self.gf('cyano.history.HistoryForeignKey')(blank=True, related_name=u'metabolites', null=True, on_delete=models.SET_NULL, to=orm['cyano.MediaComposition'])),
        ))
        db.send_create_signal(u'cyano', ['Metabolite'])

        # Adding M2M table for field biomass_composition on 'Metabolite'
        m2m_table_name = db.shorten_name(u'cyano_metabolite_biomass_composition')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('metabolite', models.ForeignKey(orm[u'cyano.metabolite'], null=False)),
            ('biomasscomposition', models.ForeignKey(orm[u'cyano.biomasscomposition'], null=False))
        ))
        db.create_unique(m2m_table_name, ['metabolite_id', 'biomasscomposition_id'])

        # Adding M2M table for field map_coordinates on 'Metabolite'
        m2m_table_name = db.shorten_name(u'cyano_metabolite_map_coordinates')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('metabolite', models.ForeignKey(orm[u'cyano.metabolite'], null=False)),
            ('metabolitemapcoordinate', models.ForeignKey(orm[u'cyano.metabolitemapcoordinate'], null=False))
        ))
        db.create_unique(m2m_table_name, ['metabolite_id', 'metabolitemapcoordinate_id'])

        # Adding model 'Note'
        db.create_table(u'cyano_note', (
            ('parent_ptr_species_component', self.gf('django.db.models.fields.related.OneToOneField')(related_name=u'child_ptr_note', unique=True, primary_key=True, to=orm['cyano.SpeciesComponent'])),
        ))
        db.send_create_signal(u'cyano', ['Note'])

        # Adding model 'Parameter'
        db.create_table(u'cyano_parameter', (
            ('parent_ptr_species_component', self.gf('django.db.models.fields.related.OneToOneField')(related_name=u'child_ptr_parameter', unique=True, primary_key=True, to=orm['cyano.SpeciesComponent'])),
            ('value', self.gf('cyano.history.HistoryForeignKey')(to=orm['cyano.EntryCharData'])),
            ('state', self.gf('cyano.history.HistoryForeignKey')(blank=True, related_name=u'parameters', null=True, on_delete=models.SET_NULL, to=orm['cyano.State'])),
            ('process', self.gf('cyano.history.HistoryForeignKey')(blank=True, related_name=u'parameters', null=True, on_delete=models.SET_NULL, to=orm['cyano.Process'])),
        ))
        db.send_create_signal(u'cyano', ['Parameter'])

        # Adding M2M table for field reactions on 'Parameter'
        m2m_table_name = db.shorten_name(u'cyano_parameter_reactions')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('parameter', models.ForeignKey(orm[u'cyano.parameter'], null=False)),
            ('reaction', models.ForeignKey(orm[u'cyano.reaction'], null=False))
        ))
        db.create_unique(m2m_table_name, ['parameter_id', 'reaction_id'])

        # Adding M2M table for field molecules on 'Parameter'
        m2m_table_name = db.shorten_name(u'cyano_parameter_molecules')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('parameter', models.ForeignKey(orm[u'cyano.parameter'], null=False)),
            ('molecule', models.ForeignKey(orm[u'cyano.molecule'], null=False))
        ))
        db.create_unique(m2m_table_name, ['parameter_id', 'molecule_id'])

        # Adding model 'Pathway'
        db.create_table(u'cyano_pathway', (
            ('parent_ptr_species_component', self.gf('django.db.models.fields.related.OneToOneField')(related_name=u'child_ptr_pathway', unique=True, primary_key=True, to=orm['cyano.SpeciesComponent'])),
        ))
        db.send_create_signal(u'cyano', ['Pathway'])

        # Adding model 'Process'
        db.create_table(u'cyano_process', (
            ('parent_ptr_species_component', self.gf('django.db.models.fields.related.OneToOneField')(related_name=u'child_ptr_process', unique=True, primary_key=True, to=orm['cyano.SpeciesComponent'])),
            ('initialization_order', self.gf('django.db.models.fields.PositiveIntegerField')()),
        ))
        db.send_create_signal(u'cyano', ['Process'])

        # Adding model 'ProteinComplex'
        db.create_table(u'cyano_proteincomplex', (
            ('parent_ptr_protein', self.gf('django.db.models.fields.related.OneToOneField')(related_name=u'child_ptr_protein_complex', unique=True, primary_key=True, to=orm['cyano.Protein'])),
            ('formation_process', self.gf('cyano.history.HistoryForeignKey')(blank=True, related_name=u'formed_complexes', null=True, on_delete=models.SET_NULL, to=orm['cyano.Process'])),
        ))
        db.send_create_signal(u'cyano', ['ProteinComplex'])

        # Adding M2M table for field biosynthesis on 'ProteinComplex'
        m2m_table_name = db.shorten_name(u'cyano_proteincomplex_biosynthesis')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('proteincomplex', models.ForeignKey(orm[u'cyano.proteincomplex'], null=False)),
            ('proteincomplexbiosythesisparticipant', models.ForeignKey(orm[u'cyano.proteincomplexbiosythesisparticipant'], null=False))
        ))
        db.create_unique(m2m_table_name, ['proteincomplex_id', 'proteincomplexbiosythesisparticipant_id'])

        # Adding M2M table for field disulfide_bonds on 'ProteinComplex'
        m2m_table_name = db.shorten_name(u'cyano_proteincomplex_disulfide_bonds')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('proteincomplex', models.ForeignKey(orm[u'cyano.proteincomplex'], null=False)),
            ('disulfidebond', models.ForeignKey(orm[u'cyano.disulfidebond'], null=False))
        ))
        db.create_unique(m2m_table_name, ['proteincomplex_id', 'disulfidebond_id'])

        # Adding model 'ProteinMonomer'
        db.create_table(u'cyano_proteinmonomer', (
            ('parent_ptr_protein', self.gf('django.db.models.fields.related.OneToOneField')(related_name=u'child_ptr_protein_monomer', unique=True, primary_key=True, to=orm['cyano.Protein'])),
            ('gene', self.gf('cyano.history.HistoryForeignKey')(related_name=u'protein_monomers', to=orm['cyano.Gene'])),
            ('is_n_terminal_methionine_cleaved', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', null=True, to=orm['cyano.EntryBooleanData'])),
            ('localization', self.gf('cyano.history.HistoryForeignKey')(related_name=u'protein_monomers', null=True, to=orm['cyano.Compartment'])),
            ('signal_sequence', self.gf('cyano.history.HistoryForeignKey')(blank=True, related_name=u'protein_monomers', null=True, on_delete=models.SET_NULL, to=orm['cyano.SignalSequence'])),
        ))
        db.send_create_signal(u'cyano', ['ProteinMonomer'])

        # Adding model 'Reaction'
        db.create_table(u'cyano_reaction', (
            ('parent_ptr_species_component', self.gf('django.db.models.fields.related.OneToOneField')(related_name=u'child_ptr_reaction', unique=True, primary_key=True, to=orm['cyano.SpeciesComponent'])),
            ('direction', self.gf('django.db.models.fields.CharField')(max_length=1)),
            ('modification', self.gf('cyano.history.HistoryForeignKey')(blank=True, related_name=u'reactions', null=True, on_delete=models.SET_NULL, to=orm['cyano.ModificationReaction'])),
            ('enzyme', self.gf('cyano.history.HistoryForeignKey')(blank=True, related_name=u'reactions', null=True, on_delete=models.SET_NULL, to=orm['cyano.EnzymeParticipant'])),
            ('is_spontaneous', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('delta_g', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('keq', self.gf('cyano.history.HistoryForeignKey')(blank=True, related_name=u'+', null=True, on_delete=models.SET_NULL, to=orm['cyano.EntryPositiveFloatData'])),
            ('kinetics_forward', self.gf('cyano.history.HistoryForeignKey')(blank=True, related_name=u'reactions_forward', null=True, on_delete=models.SET_NULL, to=orm['cyano.Kinetics'])),
            ('kinetics_backward', self.gf('cyano.history.HistoryForeignKey')(blank=True, related_name=u'reactions_backward', null=True, on_delete=models.SET_NULL, to=orm['cyano.Kinetics'])),
            ('optimal_ph', self.gf('cyano.history.HistoryForeignKey')(blank=True, related_name=u'+', null=True, on_delete=models.SET_NULL, to=orm['cyano.EntryPositiveFloatData'])),
            ('optimal_temperature', self.gf('cyano.history.HistoryForeignKey')(blank=True, related_name=u'+', null=True, on_delete=models.SET_NULL, to=orm['cyano.EntryFloatData'])),
            ('processes', self.gf('cyano.history.HistoryForeignKey')(blank=True, related_name=u'reactions', null=True, on_delete=models.SET_NULL, to=orm['cyano.Process'])),
            ('states', self.gf('cyano.history.HistoryForeignKey')(blank=True, related_name=u'reactions', null=True, on_delete=models.SET_NULL, to=orm['cyano.State'])),
        ))
        db.send_create_signal(u'cyano', ['Reaction'])

        # Adding M2M table for field stoichiometry on 'Reaction'
        m2m_table_name = db.shorten_name(u'cyano_reaction_stoichiometry')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('reaction', models.ForeignKey(orm[u'cyano.reaction'], null=False)),
            ('reactionstoichiometryparticipant', models.ForeignKey(orm[u'cyano.reactionstoichiometryparticipant'], null=False))
        ))
        db.create_unique(m2m_table_name, ['reaction_id', 'reactionstoichiometryparticipant_id'])

        # Adding M2M table for field coenzymes on 'Reaction'
        m2m_table_name = db.shorten_name(u'cyano_reaction_coenzymes')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('reaction', models.ForeignKey(orm[u'cyano.reaction'], null=False)),
            ('coenzymeparticipant', models.ForeignKey(orm[u'cyano.coenzymeparticipant'], null=False))
        ))
        db.create_unique(m2m_table_name, ['reaction_id', 'coenzymeparticipant_id'])

        # Adding M2M table for field pathways on 'Reaction'
        m2m_table_name = db.shorten_name(u'cyano_reaction_pathways')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('reaction', models.ForeignKey(orm[u'cyano.reaction'], null=False)),
            ('pathway', models.ForeignKey(orm[u'cyano.pathway'], null=False))
        ))
        db.create_unique(m2m_table_name, ['reaction_id', 'pathway_id'])

        # Adding M2M table for field map_coordinates on 'Reaction'
        m2m_table_name = db.shorten_name(u'cyano_reaction_map_coordinates')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('reaction', models.ForeignKey(orm[u'cyano.reaction'], null=False)),
            ('reactionmapcoordinate', models.ForeignKey(orm[u'cyano.reactionmapcoordinate'], null=False))
        ))
        db.create_unique(m2m_table_name, ['reaction_id', 'reactionmapcoordinate_id'])

        # Adding model 'Species'
        db.create_table(u'cyano_species', (
            ('parent_ptr_entry', self.gf('django.db.models.fields.related.OneToOneField')(related_name=u'child_ptr_species', unique=True, primary_key=True, to=orm['cyano.Entry'])),
            ('genetic_code', self.gf('django.db.models.fields.CharField')(max_length=50)),
        ))
        db.send_create_signal(u'cyano', ['Species'])

        # Adding M2M table for field cross_references on 'Species'
        m2m_table_name = db.shorten_name(u'cyano_species_cross_references')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('species', models.ForeignKey(orm[u'cyano.species'], null=False)),
            ('crossreference', models.ForeignKey(orm[u'cyano.crossreference'], null=False))
        ))
        db.create_unique(m2m_table_name, ['species_id', 'crossreference_id'])

        # Adding M2M table for field publication_references on 'Species'
        m2m_table_name = db.shorten_name(u'cyano_species_publication_references')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('species', models.ForeignKey(orm[u'cyano.species'], null=False)),
            ('publicationreference', models.ForeignKey(orm[u'cyano.publicationreference'], null=False))
        ))
        db.create_unique(m2m_table_name, ['species_id', 'publicationreference_id'])

        # Adding model 'State'
        db.create_table(u'cyano_state', (
            ('parent_ptr_species_component', self.gf('django.db.models.fields.related.OneToOneField')(related_name=u'child_ptr_state', unique=True, primary_key=True, to=orm['cyano.SpeciesComponent'])),
        ))
        db.send_create_signal(u'cyano', ['State'])

        # Adding model 'Stimulus'
        db.create_table(u'cyano_stimulus', (
            ('parent_ptr_molecule', self.gf('django.db.models.fields.related.OneToOneField')(related_name=u'child_ptr_stimulus', unique=True, primary_key=True, to=orm['cyano.Molecule'])),
            ('value', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', to=orm['cyano.EntryFloatData'])),
        ))
        db.send_create_signal(u'cyano', ['Stimulus'])

        # Adding model 'TranscriptionUnit'
        db.create_table(u'cyano_transcriptionunit', (
            ('parent_ptr_molecule', self.gf('django.db.models.fields.related.OneToOneField')(related_name=u'child_ptr_transcription_unit', unique=True, primary_key=True, to=orm['cyano.Molecule'])),
            ('promoter_35_coordinate', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('promoter_35_length', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('promoter_10_coordinate', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('promoter_10_length', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('tss_coordinate', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
        ))
        db.send_create_signal(u'cyano', ['TranscriptionUnit'])

        # Adding M2M table for field genes on 'TranscriptionUnit'
        m2m_table_name = db.shorten_name(u'cyano_transcriptionunit_genes')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('transcriptionunit', models.ForeignKey(orm[u'cyano.transcriptionunit'], null=False)),
            ('gene', models.ForeignKey(orm[u'cyano.gene'], null=False))
        ))
        db.create_unique(m2m_table_name, ['transcriptionunit_id', 'gene_id'])

        # Adding model 'TranscriptionalRegulation'
        db.create_table(u'cyano_transcriptionalregulation', (
            ('parent_ptr_species_component', self.gf('django.db.models.fields.related.OneToOneField')(related_name=u'child_ptr_transcriptional_regulation', unique=True, primary_key=True, to=orm['cyano.SpeciesComponent'])),
            ('transcription_unit', self.gf('cyano.history.HistoryForeignKey')(related_name=u'transcriptional_regulations', to=orm['cyano.TranscriptionUnit'])),
            ('transcription_factor', self.gf('cyano.history.HistoryForeignKey')(related_name=u'transcriptional_regulations', to=orm['cyano.Protein'])),
            ('binding_site', self.gf('cyano.history.HistoryForeignKey')(blank=True, related_name=u'transcriptional_regulations', null=True, on_delete=models.SET_NULL, to=orm['cyano.BindingSite'])),
            ('affinity', self.gf('cyano.history.HistoryForeignKey')(blank=True, related_name=u'+', null=True, on_delete=models.SET_NULL, to=orm['cyano.EntryPositiveFloatData'])),
            ('activity', self.gf('cyano.history.HistoryForeignKey')(related_name=u'+', to=orm['cyano.EntryPositiveFloatData'])),
        ))
        db.send_create_signal(u'cyano', ['TranscriptionalRegulation'])

        # Adding model 'Type'
        db.create_table(u'cyano_type', (
            ('parent_ptr_species_component', self.gf('django.db.models.fields.related.OneToOneField')(related_name=u'child_ptr_type', unique=True, primary_key=True, to=orm['cyano.SpeciesComponent'])),
            ('parent', self.gf('cyano.history.HistoryForeignKey')(blank=True, related_name=u'children', null=True, on_delete=models.SET_NULL, to=orm['cyano.Type'])),
        ))
        db.send_create_signal(u'cyano', ['Type'])

        # Adding model 'CrossReference'
        db.create_table(u'cyano_crossreference', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('xid', self.gf('django.db.models.fields.CharField')(max_length=255)),
            ('source', self.gf('django.db.models.fields.CharField')(max_length=20)),
        ))
        db.send_create_signal(u'cyano', ['CrossReference'])

        # Adding model 'PublicationReference'
        db.create_table(u'cyano_publicationreference', (
            ('parent_ptr_species_component', self.gf('django.db.models.fields.related.OneToOneField')(related_name=u'child_ptr_publicationreference', unique=True, primary_key=True, to=orm['cyano.SpeciesComponent'])),
            ('authors', self.gf('django.db.models.fields.TextField')(default=u'', blank=True)),
            ('editors', self.gf('django.db.models.fields.TextField')(default=u'', blank=True)),
            ('year', self.gf('django.db.models.fields.PositiveIntegerField')(null=True, blank=True)),
            ('title', self.gf('django.db.models.fields.TextField')(default=u'', blank=True)),
            ('publication', self.gf('django.db.models.fields.CharField')(default=u'', max_length=255, blank=True)),
            ('publisher', self.gf('django.db.models.fields.CharField')(default=u'', max_length=255, blank=True)),
            ('volume', self.gf('django.db.models.fields.CharField')(default=u'', max_length=255, blank=True)),
            ('issue', self.gf('django.db.models.fields.CharField')(default=u'', max_length=255, blank=True)),
            ('pages', self.gf('django.db.models.fields.CharField')(default=u'', max_length=255, blank=True)),
        ))
        db.send_create_signal(u'cyano', ['PublicationReference'])


    def backwards(self, orm):
        # Removing unique constraint on 'TableMetaManyToMany', fields ['m2m_table', 'source_table', 'target_table']
        db.delete_unique(u'cyano_tablemetamanytomany', ['m2m_table_id', 'source_table_id', 'target_table_id'])

        # Removing unique constraint on 'TableMetaColumn', fields ['table', 'column_name']
        db.delete_unique(u'cyano_tablemetacolumn', ['table_id', 'column_name'])

        # Removing unique constraint on 'GroupPermission', fields ['entry', 'group']
        db.delete_unique(u'cyano_grouppermission', ['entry_id', 'group_id'])

        # Removing unique constraint on 'UserPermission', fields ['entry', 'user']
        db.delete_unique(u'cyano_userpermission', ['entry_id', 'user_id'])

        # Deleting model 'Permission'
        db.delete_table(u'cyano_permission')

        # Deleting model 'UserPermission'
        db.delete_table(u'cyano_userpermission')

        # Removing M2M table for field allow on 'UserPermission'
        db.delete_table(db.shorten_name(u'cyano_userpermission_allow'))

        # Removing M2M table for field deny on 'UserPermission'
        db.delete_table(db.shorten_name(u'cyano_userpermission_deny'))

        # Deleting model 'GroupPermission'
        db.delete_table(u'cyano_grouppermission')

        # Removing M2M table for field allow on 'GroupPermission'
        db.delete_table(db.shorten_name(u'cyano_grouppermission_allow'))

        # Removing M2M table for field deny on 'GroupPermission'
        db.delete_table(db.shorten_name(u'cyano_grouppermission_deny'))

        # Deleting model 'GroupProfile'
        db.delete_table(u'cyano_groupprofile')

        # Deleting model 'UserProfile'
        db.delete_table(u'cyano_userprofile')

        # Deleting model 'TableMeta'
        db.delete_table(u'cyano_tablemeta')

        # Deleting model 'TableMetaColumn'
        db.delete_table(u'cyano_tablemetacolumn')

        # Deleting model 'TableMetaManyToMany'
        db.delete_table(u'cyano_tablemetamanytomany')

        # Deleting model 'RevisionDetail'
        db.delete_table(u'cyano_revisiondetail')

        # Deleting model 'Revision'
        db.delete_table(u'cyano_revision')

        # Deleting model 'RevisionManyToMany'
        db.delete_table(u'cyano_revisionmanytomany')

        # Deleting model 'Evidence'
        db.delete_table(u'cyano_evidence')

        # Removing M2M table for field species_component on 'Evidence'
        db.delete_table(db.shorten_name(u'cyano_evidence_species_component'))

        # Deleting model 'EntryBooleanData'
        db.delete_table(u'cyano_entrybooleandata')

        # Removing M2M table for field evidence on 'EntryBooleanData'
        db.delete_table(db.shorten_name(u'cyano_entrybooleandata_evidence'))

        # Deleting model 'EntryCharData'
        db.delete_table(u'cyano_entrychardata')

        # Removing M2M table for field evidence on 'EntryCharData'
        db.delete_table(db.shorten_name(u'cyano_entrychardata_evidence'))

        # Deleting model 'EntryFloatData'
        db.delete_table(u'cyano_entryfloatdata')

        # Removing M2M table for field evidence on 'EntryFloatData'
        db.delete_table(db.shorten_name(u'cyano_entryfloatdata_evidence'))

        # Deleting model 'EntryPositiveFloatData'
        db.delete_table(u'cyano_entrypositivefloatdata')

        # Removing M2M table for field evidence on 'EntryPositiveFloatData'
        db.delete_table(db.shorten_name(u'cyano_entrypositivefloatdata_evidence'))

        # Deleting model 'EntryTextData'
        db.delete_table(u'cyano_entrytextdata')

        # Removing M2M table for field evidence on 'EntryTextData'
        db.delete_table(db.shorten_name(u'cyano_entrytextdata_evidence'))

        # Deleting model 'BindingSite'
        db.delete_table(u'cyano_bindingsite')

        # Removing M2M table for field evidence on 'BindingSite'
        db.delete_table(db.shorten_name(u'cyano_bindingsite_evidence'))

        # Deleting model 'BiomassComposition'
        db.delete_table(u'cyano_biomasscomposition')

        # Removing M2M table for field evidence on 'BiomassComposition'
        db.delete_table(db.shorten_name(u'cyano_biomasscomposition_evidence'))

        # Deleting model 'Codon'
        db.delete_table(u'cyano_codon')

        # Removing M2M table for field evidence on 'Codon'
        db.delete_table(db.shorten_name(u'cyano_codon_evidence'))

        # Deleting model 'CoenzymeParticipant'
        db.delete_table(u'cyano_coenzymeparticipant')

        # Removing M2M table for field evidence on 'CoenzymeParticipant'
        db.delete_table(db.shorten_name(u'cyano_coenzymeparticipant_evidence'))

        # Deleting model 'DisulfideBond'
        db.delete_table(u'cyano_disulfidebond')

        # Removing M2M table for field evidence on 'DisulfideBond'
        db.delete_table(db.shorten_name(u'cyano_disulfidebond_evidence'))

        # Deleting model 'DNAFootprint'
        db.delete_table(u'cyano_dnafootprint')

        # Removing M2M table for field evidence on 'DNAFootprint'
        db.delete_table(db.shorten_name(u'cyano_dnafootprint_evidence'))

        # Deleting model 'EnzymeParticipant'
        db.delete_table(u'cyano_enzymeparticipant')

        # Removing M2M table for field evidence on 'EnzymeParticipant'
        db.delete_table(db.shorten_name(u'cyano_enzymeparticipant_evidence'))

        # Deleting model 'Homolog'
        db.delete_table(u'cyano_homolog')

        # Removing M2M table for field evidence on 'Homolog'
        db.delete_table(db.shorten_name(u'cyano_homolog_evidence'))

        # Deleting model 'Kinetics'
        db.delete_table(u'cyano_kinetics')

        # Removing M2M table for field evidence on 'Kinetics'
        db.delete_table(db.shorten_name(u'cyano_kinetics_evidence'))

        # Deleting model 'MediaComposition'
        db.delete_table(u'cyano_mediacomposition')

        # Removing M2M table for field evidence on 'MediaComposition'
        db.delete_table(db.shorten_name(u'cyano_mediacomposition_evidence'))

        # Deleting model 'MetaboliteMapCoordinate'
        db.delete_table(u'cyano_metabolitemapcoordinate')

        # Deleting model 'ModificationReaction'
        db.delete_table(u'cyano_modificationreaction')

        # Removing M2M table for field evidence on 'ModificationReaction'
        db.delete_table(db.shorten_name(u'cyano_modificationreaction_evidence'))

        # Deleting model 'ProstheticGroupParticipant'
        db.delete_table(u'cyano_prostheticgroupparticipant')

        # Removing M2M table for field evidence on 'ProstheticGroupParticipant'
        db.delete_table(db.shorten_name(u'cyano_prostheticgroupparticipant_evidence'))

        # Deleting model 'ProteinComplexBiosythesisParticipant'
        db.delete_table(u'cyano_proteincomplexbiosythesisparticipant')

        # Removing M2M table for field evidence on 'ProteinComplexBiosythesisParticipant'
        db.delete_table(db.shorten_name(u'cyano_proteincomplexbiosythesisparticipant_evidence'))

        # Deleting model 'ReactionMapCoordinate'
        db.delete_table(u'cyano_reactionmapcoordinate')

        # Deleting model 'ReactionStoichiometryParticipant'
        db.delete_table(u'cyano_reactionstoichiometryparticipant')

        # Removing M2M table for field evidence on 'ReactionStoichiometryParticipant'
        db.delete_table(db.shorten_name(u'cyano_reactionstoichiometryparticipant_evidence'))

        # Deleting model 'SignalSequence'
        db.delete_table(u'cyano_signalsequence')

        # Removing M2M table for field evidence on 'SignalSequence'
        db.delete_table(db.shorten_name(u'cyano_signalsequence_evidence'))

        # Deleting model 'Synonym'
        db.delete_table(u'cyano_synonym')

        # Deleting model 'Entry'
        db.delete_table(u'cyano_entry')

        # Removing M2M table for field synonyms on 'Entry'
        db.delete_table(db.shorten_name(u'cyano_entry_synonyms'))

        # Deleting model 'CrossReferenceMeta'
        db.delete_table(u'cyano_crossreferencemeta')

        # Deleting model 'SpeciesComponent'
        db.delete_table(u'cyano_speciescomponent')

        # Removing M2M table for field species on 'SpeciesComponent'
        db.delete_table(db.shorten_name(u'cyano_speciescomponent_species'))

        # Removing M2M table for field type on 'SpeciesComponent'
        db.delete_table(db.shorten_name(u'cyano_speciescomponent_type'))

        # Removing M2M table for field cross_references on 'SpeciesComponent'
        db.delete_table(db.shorten_name(u'cyano_speciescomponent_cross_references'))

        # Removing M2M table for field publication_references on 'SpeciesComponent'
        db.delete_table(db.shorten_name(u'cyano_speciescomponent_publication_references'))

        # Deleting model 'Molecule'
        db.delete_table(u'cyano_molecule')

        # Deleting model 'Protein'
        db.delete_table(u'cyano_protein')

        # Removing M2M table for field prosthetic_groups on 'Protein'
        db.delete_table(db.shorten_name(u'cyano_protein_prosthetic_groups'))

        # Removing M2M table for field chaperones on 'Protein'
        db.delete_table(db.shorten_name(u'cyano_protein_chaperones'))

        # Deleting model 'DNA'
        db.delete_table(u'cyano_dna')

        # Deleting model 'Chromosome'
        db.delete_table(u'cyano_chromosome')

        # Deleting model 'Plasmid'
        db.delete_table(u'cyano_plasmid')

        # Deleting model 'ChromosomeFeature'
        db.delete_table(u'cyano_chromosomefeature')

        # Deleting model 'Compartment'
        db.delete_table(u'cyano_compartment')

        # Deleting model 'Gene'
        db.delete_table(u'cyano_gene')

        # Removing M2M table for field codons on 'Gene'
        db.delete_table(db.shorten_name(u'cyano_gene_codons'))

        # Removing M2M table for field homologs on 'Gene'
        db.delete_table(db.shorten_name(u'cyano_gene_homologs'))

        # Deleting model 'Metabolite'
        db.delete_table(u'cyano_metabolite')

        # Removing M2M table for field biomass_composition on 'Metabolite'
        db.delete_table(db.shorten_name(u'cyano_metabolite_biomass_composition'))

        # Removing M2M table for field map_coordinates on 'Metabolite'
        db.delete_table(db.shorten_name(u'cyano_metabolite_map_coordinates'))

        # Deleting model 'Note'
        db.delete_table(u'cyano_note')

        # Deleting model 'Parameter'
        db.delete_table(u'cyano_parameter')

        # Removing M2M table for field reactions on 'Parameter'
        db.delete_table(db.shorten_name(u'cyano_parameter_reactions'))

        # Removing M2M table for field molecules on 'Parameter'
        db.delete_table(db.shorten_name(u'cyano_parameter_molecules'))

        # Deleting model 'Pathway'
        db.delete_table(u'cyano_pathway')

        # Deleting model 'Process'
        db.delete_table(u'cyano_process')

        # Deleting model 'ProteinComplex'
        db.delete_table(u'cyano_proteincomplex')

        # Removing M2M table for field biosynthesis on 'ProteinComplex'
        db.delete_table(db.shorten_name(u'cyano_proteincomplex_biosynthesis'))

        # Removing M2M table for field disulfide_bonds on 'ProteinComplex'
        db.delete_table(db.shorten_name(u'cyano_proteincomplex_disulfide_bonds'))

        # Deleting model 'ProteinMonomer'
        db.delete_table(u'cyano_proteinmonomer')

        # Deleting model 'Reaction'
        db.delete_table(u'cyano_reaction')

        # Removing M2M table for field stoichiometry on 'Reaction'
        db.delete_table(db.shorten_name(u'cyano_reaction_stoichiometry'))

        # Removing M2M table for field coenzymes on 'Reaction'
        db.delete_table(db.shorten_name(u'cyano_reaction_coenzymes'))

        # Removing M2M table for field pathways on 'Reaction'
        db.delete_table(db.shorten_name(u'cyano_reaction_pathways'))

        # Removing M2M table for field map_coordinates on 'Reaction'
        db.delete_table(db.shorten_name(u'cyano_reaction_map_coordinates'))

        # Deleting model 'Species'
        db.delete_table(u'cyano_species')

        # Removing M2M table for field cross_references on 'Species'
        db.delete_table(db.shorten_name(u'cyano_species_cross_references'))

        # Removing M2M table for field publication_references on 'Species'
        db.delete_table(db.shorten_name(u'cyano_species_publication_references'))

        # Deleting model 'State'
        db.delete_table(u'cyano_state')

        # Deleting model 'Stimulus'
        db.delete_table(u'cyano_stimulus')

        # Deleting model 'TranscriptionUnit'
        db.delete_table(u'cyano_transcriptionunit')

        # Removing M2M table for field genes on 'TranscriptionUnit'
        db.delete_table(db.shorten_name(u'cyano_transcriptionunit_genes'))

        # Deleting model 'TranscriptionalRegulation'
        db.delete_table(u'cyano_transcriptionalregulation')

        # Deleting model 'Type'
        db.delete_table(u'cyano_type')

        # Deleting model 'CrossReference'
        db.delete_table(u'cyano_crossreference')

        # Deleting model 'PublicationReference'
        db.delete_table(u'cyano_publicationreference')


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
            'groups': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['auth.Group']", 'symmetrical': 'False', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'is_active': ('django.db.models.fields.BooleanField', [], {'default': 'True'}),
            'is_staff': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'is_superuser': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'last_login': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'last_name': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'password': ('django.db.models.fields.CharField', [], {'max_length': '128'}),
            'user_permissions': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['auth.Permission']", 'symmetrical': 'False', 'blank': 'True'}),
            'username': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '30'})
        },
        u'contenttypes.contenttype': {
            'Meta': {'ordering': "('name',)", 'unique_together': "(('app_label', 'model'),)", 'object_name': 'ContentType', 'db_table': "'django_content_type'"},
            'app_label': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'model': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100'})
        },
        u'cyano.bindingsite': {
            'Meta': {'ordering': "[u'coordinate', u'length']", 'object_name': 'BindingSite'},
            'coordinate': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'detail': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.RevisionDetail']"}),
            'direction': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'length': ('django.db.models.fields.PositiveIntegerField', [], {})
        },
        u'cyano.biomasscomposition': {
            'Meta': {'ordering': "[u'-concentration']", 'object_name': 'BiomassComposition'},
            'compartment': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'biomass_compositions'", 'to': u"orm['cyano.Compartment']"}),
            'concentration': ('django.db.models.fields.FloatField', [], {}),
            'detail': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.RevisionDetail']"}),
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'})
        },
        u'cyano.chromosome': {
            'Meta': {'object_name': 'Chromosome', '_ormbases': [u'cyano.DNA']},
            'parent_ptr_dna': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "u'child_ptr_chromosome'", 'unique': 'True', 'primary_key': 'True', 'to': u"orm['cyano.DNA']"})
        },
        u'cyano.chromosomefeature': {
            'Meta': {'object_name': 'ChromosomeFeature', '_ormbases': [u'cyano.SpeciesComponent']},
            'chromosome': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'features'", 'to': u"orm['cyano.DNA']"}),
            'coordinate': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'direction': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'length': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'parent_ptr_species_component': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "u'child_ptr_chromosome_feature'", 'unique': 'True', 'primary_key': 'True', 'to': u"orm['cyano.SpeciesComponent']"})
        },
        u'cyano.codon': {
            'Meta': {'ordering': "[u'sequence']", 'object_name': 'Codon'},
            'detail': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.RevisionDetail']"}),
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'sequence': ('django.db.models.fields.CharField', [], {'max_length': '3'})
        },
        u'cyano.coenzymeparticipant': {
            'Meta': {'object_name': 'CoenzymeParticipant'},
            'coefficient': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'compartment': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.Compartment']"}),
            'detail': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.RevisionDetail']"}),
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'metabolite': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'coenzyme_participants'", 'to': u"orm['cyano.Metabolite']"})
        },
        u'cyano.compartment': {
            'Meta': {'object_name': 'Compartment', '_ormbases': [u'cyano.SpeciesComponent']},
            'parent_ptr_species_component': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "u'child_ptr_compartment'", 'unique': 'True', 'primary_key': 'True', 'to': u"orm['cyano.SpeciesComponent']"})
        },
        u'cyano.crossreference': {
            'Meta': {'ordering': "[u'xid']", 'object_name': 'CrossReference'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'source': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'xid': ('django.db.models.fields.CharField', [], {'max_length': '255'})
        },
        u'cyano.crossreferencemeta': {
            'Meta': {'object_name': 'CrossReferenceMeta'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '255'})
        },
        u'cyano.disulfidebond': {
            'Meta': {'object_name': 'DisulfideBond'},
            'detail': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.RevisionDetail']"}),
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'protein_monomer': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'disulfide_bonds'", 'to': u"orm['cyano.ProteinMonomer']"}),
            'residue_1': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'residue_2': ('django.db.models.fields.PositiveIntegerField', [], {})
        },
        u'cyano.dna': {
            'Meta': {'object_name': 'DNA', '_ormbases': [u'cyano.Molecule']},
            'length': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'parent_ptr_molecule': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "u'child_ptr_dna'", 'unique': 'True', 'primary_key': 'True', 'to': u"orm['cyano.Molecule']"}),
            'sequence': ('django.db.models.fields.TextField', [], {'default': "u''", 'blank': 'True'})
        },
        u'cyano.dnafootprint': {
            'Meta': {'object_name': 'DNAFootprint'},
            'binding': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '10', 'blank': 'True'}),
            'detail': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.RevisionDetail']"}),
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'length': ('django.db.models.fields.PositiveIntegerField', [], {'null': 'True', 'blank': 'True'}),
            'region': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '10', 'blank': 'True'})
        },
        u'cyano.entry': {
            'Meta': {'ordering': "[u'wid']", 'object_name': 'Entry'},
            'comments': ('django.db.models.fields.TextField', [], {'default': "u''", 'blank': 'True'}),
            'created_detail': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'entry_created_detail'", 'to': u"orm['cyano.RevisionDetail']"}),
            'detail': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'entry_detail'", 'to': u"orm['cyano.RevisionDetail']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'model_type': ('cyano.history.HistoryForeignKey', [], {'to': u"orm['cyano.TableMeta']"}),
            'name': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255', 'blank': 'True'}),
            'synonyms': ('cyano.history.HistoryManyToManyField', [], {'blank': 'True', 'related_name': "u'entry'", 'null': 'True', 'symmetrical': 'False', 'to': u"orm['cyano.Synonym']"}),
            'wid': ('django.db.models.fields.SlugField', [], {'max_length': '150'})
        },
        u'cyano.entrybooleandata': {
            'Meta': {'ordering': "[u'value']", 'object_name': 'EntryBooleanData'},
            'detail': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.RevisionDetail']"}),
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'value': ('django.db.models.fields.BooleanField', [], {'default': 'False'})
        },
        u'cyano.entrychardata': {
            'Meta': {'ordering': "[u'value', u'units']", 'object_name': 'EntryCharData'},
            'detail': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.RevisionDetail']"}),
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'units': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255', 'blank': 'True'}),
            'value': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255', 'blank': 'True'})
        },
        u'cyano.entryfloatdata': {
            'Meta': {'ordering': "[u'value', u'units']", 'object_name': 'EntryFloatData'},
            'detail': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.RevisionDetail']"}),
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'units': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255', 'blank': 'True'}),
            'value': ('django.db.models.fields.FloatField', [], {})
        },
        u'cyano.entrypositivefloatdata': {
            'Meta': {'ordering': "[u'value', u'units']", 'object_name': 'EntryPositiveFloatData'},
            'detail': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.RevisionDetail']"}),
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'units': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255', 'blank': 'True'}),
            'value': ('django.db.models.fields.FloatField', [], {})
        },
        u'cyano.entrytextdata': {
            'Meta': {'ordering': "[u'value', u'units']", 'object_name': 'EntryTextData'},
            'detail': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.RevisionDetail']"}),
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'units': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255', 'blank': 'True'}),
            'value': ('django.db.models.fields.TextField', [], {'default': "u''", 'blank': 'True'})
        },
        u'cyano.enzymeparticipant': {
            'Meta': {'object_name': 'EnzymeParticipant'},
            'compartment': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.Compartment']"}),
            'detail': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.RevisionDetail']"}),
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'protein': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'enzyme_participants'", 'to': u"orm['cyano.Protein']"})
        },
        u'cyano.evidence': {
            'Meta': {'ordering': "[u'value', u'units']", 'object_name': 'Evidence'},
            'comments': ('django.db.models.fields.TextField', [], {'default': "u''", 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'is_experimentally_constrained': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'media': ('django.db.models.fields.CharField', [], {'max_length': '255', 'null': 'True', 'blank': 'True'}),
            'pH': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'references': ('cyano.history.HistoryForeignKey', [], {'blank': 'True', 'related_name': "u'evidence'", 'null': 'True', 'to': u"orm['cyano.PublicationReference']"}),
            'species': ('django.db.models.fields.CharField', [], {'max_length': '255', 'null': 'True', 'blank': 'True'}),
            'species_component': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'related_name': "u'+'", 'blank': 'True', 'to': u"orm['cyano.SpeciesComponent']"}),
            'temperature': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'units': ('django.db.models.fields.CharField', [], {'max_length': '255', 'null': 'True', 'blank': 'True'}),
            'value': ('django.db.models.fields.TextField', [], {'default': "u''", 'blank': 'True'})
        },
        u'cyano.gene': {
            'Meta': {'object_name': 'Gene', '_ormbases': [u'cyano.Molecule']},
            'amino_acid': ('cyano.history.HistoryForeignKey', [], {'blank': 'True', 'related_name': "u'genes'", 'null': 'True', 'on_delete': 'models.SET_NULL', 'to': u"orm['cyano.Metabolite']"}),
            'chromosome': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'genes'", 'to': u"orm['cyano.DNA']"}),
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
            'detail': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.RevisionDetail']"}),
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'species': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'xid': ('django.db.models.fields.CharField', [], {'max_length': '255'})
        },
        u'cyano.kinetics': {
            'Meta': {'object_name': 'Kinetics'},
            'detail': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.RevisionDetail']"}),
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'km': ('django.db.models.fields.CharField', [], {'max_length': '255', 'blank': 'True'}),
            'rate_law': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255', 'blank': 'True'}),
            'vmax': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'vmax_unit': ('django.db.models.fields.CharField', [], {'max_length': '255', 'blank': 'True'})
        },
        u'cyano.mediacomposition': {
            'Meta': {'ordering': "[u'-concentration']", 'object_name': 'MediaComposition'},
            'concentration': ('django.db.models.fields.FloatField', [], {}),
            'detail': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.RevisionDetail']"}),
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'is_diffused': ('django.db.models.fields.BooleanField', [], {'default': 'False'})
        },
        u'cyano.metabolite': {
            'Meta': {'object_name': 'Metabolite', '_ormbases': [u'cyano.Molecule']},
            'biomass_composition': ('cyano.history.HistoryManyToManyField', [], {'blank': 'True', 'related_name': "u'metabolites'", 'null': 'True', 'symmetrical': 'False', 'to': u"orm['cyano.BiomassComposition']"}),
            'charge': ('django.db.models.fields.IntegerField', [], {}),
            'deltag_formation': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'empirical_formula': ('django.db.models.fields.TextField', [], {}),
            'is_hydrophobic': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
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
            'detail': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.RevisionDetail']"}),
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
        u'cyano.permission': {
            'Meta': {'object_name': 'Permission'},
            'description': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255'})
        },
        u'cyano.plasmid': {
            'Meta': {'object_name': 'Plasmid', '_ormbases': [u'cyano.DNA']},
            'parent_ptr_dna': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "u'child_ptr_plasmid'", 'unique': 'True', 'primary_key': 'True', 'to': u"orm['cyano.DNA']"})
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
            'detail': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.RevisionDetail']"}),
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
            'detail': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.RevisionDetail']"}),
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
            'is_spontaneous': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
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
            'detail': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.RevisionDetail']"}),
            'evidence': ('cyano.history.HistoryManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['cyano.Evidence']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'molecule': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'reaction_stoichiometry_participants'", 'to': u"orm['cyano.Molecule']"})
        },
        u'cyano.revision': {
            'Meta': {'object_name': 'Revision'},
            'action': ('django.db.models.fields.CharField', [], {'max_length': '1'}),
            'column': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.TableMetaColumn']"}),
            'current': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'revisions'", 'to': u"orm['cyano.Entry']"}),
            'detail': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'revisions'", 'null': 'True', 'to': u"orm['cyano.RevisionDetail']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'new_value': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'})
        },
        u'cyano.revisiondetail': {
            'Meta': {'object_name': 'RevisionDetail'},
            'date': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'reason': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255', 'blank': 'True'}),
            'user': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.UserProfile']"})
        },
        u'cyano.revisionmanytomany': {
            'Meta': {'object_name': 'RevisionManyToMany'},
            'action': ('django.db.models.fields.CharField', [], {'max_length': '1'}),
            'current': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'revisions_m2m'", 'to': u"orm['cyano.Entry']"}),
            'detail': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'revisions_m2m'", 'null': 'True', 'to': u"orm['cyano.RevisionDetail']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'new_value': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'table': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.TableMetaManyToMany']"})
        },
        u'cyano.signalsequence': {
            'Meta': {'ordering': "[u'type', u'location', u'length']", 'object_name': 'SignalSequence'},
            'detail': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'to': u"orm['cyano.RevisionDetail']"}),
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
        u'cyano.tablemetacolumn': {
            'Meta': {'unique_together': "((u'table', u'column_name'),)", 'object_name': 'TableMetaColumn'},
            'column_id': ('django.db.models.fields.IntegerField', [], {}),
            'column_name': ('django.db.models.fields.CharField', [], {'max_length': '255'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'table': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'columns'", 'to': u"orm['cyano.TableMeta']"})
        },
        u'cyano.tablemetamanytomany': {
            'Meta': {'unique_together': "((u'm2m_table', u'source_table', u'target_table'),)", 'object_name': 'TableMetaManyToMany'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'm2m_table': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'+'", 'unique': 'True', 'to': u"orm['cyano.TableMeta']"}),
            'source_table': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'm2ms_source'", 'to': u"orm['cyano.TableMeta']"}),
            'target_table': ('cyano.history.HistoryForeignKey', [], {'related_name': "u'm2ms_target'", 'to': u"orm['cyano.TableMeta']"})
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
            'parent': ('cyano.history.HistoryForeignKey', [], {'blank': 'True', 'related_name': "u'children'", 'null': 'True', 'on_delete': 'models.SET_NULL', 'to': u"orm['cyano.Type']"}),
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