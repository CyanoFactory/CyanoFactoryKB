# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cyano', '0001_initial'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='globalpermission',
            options={'permissions': (('access_species', 'Can access any species'), ('create_mutant', 'Can create new species or mutants'), ('access_sbgn', 'Can access SBGN map'))},
        ),
        migrations.RemoveField(
            model_name='massspectrometryprotein',
            name='protein_ptr',
        ),
        migrations.RemoveField(
            model_name='peptide',
            name='protein_ptr',
        ),
        migrations.RenameField(
            model_name='massspectrometryprotein',
            old_name='parent_ptr_species_component_id',
            new_name='parent_ptr_protein'
        ),
        migrations.RenameField(
            model_name='peptide',
            old_name='parent_ptr_species_component_id',
            new_name='parent_ptr_protein'
        ),
        migrations.AlterField(
            model_name='massspectrometryprotein',
            name='parent_ptr_protein',
            field=models.OneToOneField(parent_link=True, related_name='child_ptr_ms_protein', primary_key=True, serialize=False, to='cyano.Protein', verbose_name='Protein'),
        ),
        migrations.AlterField(
            model_name='peptide',
            name='parent_ptr_protein',
            field=models.OneToOneField(parent_link=True, related_name='child_ptr_peptide', primary_key=True, serialize=False, to='cyano.Protein', verbose_name='Species component'),
        ),
    ]
