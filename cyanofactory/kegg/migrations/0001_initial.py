# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cyano', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='EcNumber',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=255)),
            ],
        ),
        migrations.CreateModel(
            name='GlobalPermission',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
            ],
            options={
                'permissions': (('access_kegg', 'Can access Kegg maps'),),
            },
        ),
        migrations.CreateModel(
            name='Map',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=255)),
                ('title', models.CharField(max_length=255, blank=True)),
                ('overview', models.BooleanField(default=False, verbose_name=b'Is an overview pathway map')),
                ('ec_numbers', models.ManyToManyField(to='kegg.EcNumber', verbose_name=b'Ec numbers belonging to this map')),
            ],
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
    ]
