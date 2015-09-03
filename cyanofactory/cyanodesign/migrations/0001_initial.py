# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import datetime
import django_extensions.db.fields.json


class Migration(migrations.Migration):

    dependencies = [
        ('cyano', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='DesignModel',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=255, verbose_name=b'Name')),
                ('filename', models.CharField(max_length=255, verbose_name=b'Filename')),
                ('content', models.TextField(verbose_name=b'BioOpt file content (deprecated)')),
                ('user', models.ForeignKey(related_name='+', editable=False, to='cyano.UserProfile', verbose_name=b'Saved by')),
            ],
        ),
        migrations.CreateModel(
            name='DesignTemplate',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=255, verbose_name=b'Name')),
                ('description', models.TextField(verbose_name=b'Description', blank=True)),
                ('filename', models.CharField(max_length=255, verbose_name=b'Filename')),
                ('content', django_extensions.db.fields.json.JSONField(default=b'{}', verbose_name=b'BioOpt file content in json format')),
            ],
        ),
        migrations.CreateModel(
            name='GlobalPermission',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
            ],
            options={
                'permissions': (('access_cyanodesign', 'Can access CyanoDesign'),),
            },
        ),
        migrations.CreateModel(
            name='Revision',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('content', django_extensions.db.fields.json.JSONField(default=b'{}', verbose_name=b'BioOpt file content in json format')),
                ('date', models.DateTimeField(default=datetime.datetime.now, verbose_name=b'Modification date')),
                ('changes', django_extensions.db.fields.json.JSONField(default=b'{}', verbose_name=b'Summary of changes')),
                ('reason', models.TextField(default=b'', verbose_name=b'Description of changes', blank=True)),
                ('model', models.ForeignKey(related_name='revisions', verbose_name=b'Model', to='cyanodesign.DesignModel')),
            ],
            options={
                'ordering': ['-date'],
            },
        ),
    ]
