# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cyano', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='BioMolecule',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('x', models.IntegerField()),
                ('y', models.IntegerField()),
                ('w', models.IntegerField()),
                ('h', models.IntegerField()),
                ('title', models.CharField(default=b'', max_length=255, verbose_name=b'title', blank=True)),
            ],
        ),
        migrations.CreateModel(
            name='Color',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(default=b'', max_length=255, verbose_name=b'color', blank=True)),
            ],
        ),
        migrations.CreateModel(
            name='GlobalPermission',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
            ],
            options={
                'permissions': (('access_boehringer', 'Can access Boehringer map'),),
            },
        ),
        migrations.CreateModel(
            name='Query',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=255, verbose_name=b'Name of the query')),
                ('query', models.TextField(verbose_name=b'Query text')),
                ('user', models.ForeignKey(related_name='+', editable=False, to='cyano.UserProfile', verbose_name=b'Saved by')),
            ],
        ),
        migrations.CreateModel(
            name='Enzyme',
            fields=[
                ('biomolecule_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='boehringer.BioMolecule')),
                ('ec', models.CharField(default=b'', max_length=255, verbose_name=b'ec-number', blank=True)),
            ],
            bases=('boehringer.biomolecule',),
        ),
        migrations.CreateModel(
            name='Metabolite',
            fields=[
                ('biomolecule_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='boehringer.BioMolecule')),
            ],
            bases=('boehringer.biomolecule',),
        ),
        migrations.AddField(
            model_name='biomolecule',
            name='color',
            field=models.ForeignKey(to='boehringer.Color'),
        ),
    ]
