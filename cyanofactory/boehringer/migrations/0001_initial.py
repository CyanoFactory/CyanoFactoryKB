# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'Color'
        db.create_table(u'boehringer_color', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(default='', max_length=255, blank=True)),
        ))
        db.send_create_signal(u'boehringer', ['Color'])

        # Adding model 'BioMolecule'
        db.create_table(u'boehringer_biomolecule', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('x', self.gf('django.db.models.fields.IntegerField')()),
            ('y', self.gf('django.db.models.fields.IntegerField')()),
            ('w', self.gf('django.db.models.fields.IntegerField')()),
            ('h', self.gf('django.db.models.fields.IntegerField')()),
            ('color', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['boehringer.Color'])),
            ('title', self.gf('django.db.models.fields.CharField')(default='', max_length=255, blank=True)),
        ))
        db.send_create_signal(u'boehringer', ['BioMolecule'])

        # Adding model 'Enzyme'
        db.create_table(u'boehringer_enzyme', (
            (u'biomolecule_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['boehringer.BioMolecule'], unique=True, primary_key=True)),
            ('ec', self.gf('django.db.models.fields.CharField')(default='', max_length=255, blank=True)),
        ))
        db.send_create_signal(u'boehringer', ['Enzyme'])

        # Adding model 'Metabolite'
        db.create_table(u'boehringer_metabolite', (
            (u'biomolecule_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['boehringer.BioMolecule'], unique=True, primary_key=True)),
        ))
        db.send_create_signal(u'boehringer', ['Metabolite'])

        # Adding model 'Query'
        db.create_table(u'boehringer_query', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('user', self.gf('django.db.models.fields.related.ForeignKey')(related_name='+', to=orm['cyano.UserProfile'])),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=255)),
            ('query', self.gf('django.db.models.fields.TextField')()),
        ))
        db.send_create_signal(u'boehringer', ['Query'])


    def backwards(self, orm):
        # Deleting model 'Color'
        db.delete_table(u'boehringer_color')

        # Deleting model 'BioMolecule'
        db.delete_table(u'boehringer_biomolecule')

        # Deleting model 'Enzyme'
        db.delete_table(u'boehringer_enzyme')

        # Deleting model 'Metabolite'
        db.delete_table(u'boehringer_metabolite')

        # Deleting model 'Query'
        db.delete_table(u'boehringer_query')


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
        u'boehringer.biomolecule': {
            'Meta': {'object_name': 'BioMolecule'},
            'color': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['boehringer.Color']"}),
            'h': ('django.db.models.fields.IntegerField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'title': ('django.db.models.fields.CharField', [], {'default': "''", 'max_length': '255', 'blank': 'True'}),
            'w': ('django.db.models.fields.IntegerField', [], {}),
            'x': ('django.db.models.fields.IntegerField', [], {}),
            'y': ('django.db.models.fields.IntegerField', [], {})
        },
        u'boehringer.color': {
            'Meta': {'object_name': 'Color'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'default': "''", 'max_length': '255', 'blank': 'True'})
        },
        u'boehringer.enzyme': {
            'Meta': {'object_name': 'Enzyme', '_ormbases': [u'boehringer.BioMolecule']},
            u'biomolecule_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': u"orm['boehringer.BioMolecule']", 'unique': 'True', 'primary_key': 'True'}),
            'ec': ('django.db.models.fields.CharField', [], {'default': "''", 'max_length': '255', 'blank': 'True'})
        },
        u'boehringer.metabolite': {
            'Meta': {'object_name': 'Metabolite', '_ormbases': [u'boehringer.BioMolecule']},
            u'biomolecule_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': u"orm['boehringer.BioMolecule']", 'unique': 'True', 'primary_key': 'True'})
        },
        u'boehringer.query': {
            'Meta': {'object_name': 'Query'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '255'}),
            'query': ('django.db.models.fields.TextField', [], {}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'+'", 'to': u"orm['cyano.UserProfile']"})
        },
        u'contenttypes.contenttype': {
            'Meta': {'ordering': "('name',)", 'unique_together': "(('app_label', 'model'),)", 'object_name': 'ContentType', 'db_table': "'django_content_type'"},
            'app_label': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'model': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100'})
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

    complete_apps = ['boehringer']