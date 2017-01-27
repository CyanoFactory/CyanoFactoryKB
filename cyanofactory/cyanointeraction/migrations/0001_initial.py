# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models
from django.contrib.postgres import fields

class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Abstracts',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('abstract_id', models.CharField(max_length=250)),
                ('publication_date', models.SmallIntegerField(null=True, blank=True)),
                ('publication_source', models.TextField(blank=True)),
                ('linkout_url', models.TextField(blank=True)),
                ('title', models.TextField()),
                ('body', models.TextField()),
            ],
            options={
                'db_table': 'abstracts',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='Actions',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('item_id_a', models.IntegerField()),
                ('item_id_b', models.IntegerField()),
                ('mode', models.CharField(max_length=255)),
                ('action', models.CharField(max_length=255)),
                ('a_is_acting', models.BooleanField()),
                ('score', models.SmallIntegerField()),
            ],
            options={
                'db_table': 'actions',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='ActionsSets',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('item_id_a', models.IntegerField()),
                ('item_id_b', models.IntegerField()),
                ('mode', models.CharField(max_length=255)),
                ('sources', models.TextField()),
            ],
            options={
                'db_table': 'actions_sets',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='BestCombinedScoresOrthgroups',
            fields=[
                ('orthgroup_id', models.IntegerField(serialize=False, primary_key=True)),
                ('best_score', models.IntegerField()),
            ],
            options={
                'db_table': 'best_combined_scores_orthgroups',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='BestCombinedScoresProteins',
            fields=[
                ('protein_id', models.IntegerField(serialize=False, primary_key=True)),
                ('best_score', models.IntegerField()),
            ],
            options={
                'db_table': 'best_combined_scores_proteins',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='ChemicalAliases',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('chemical', models.CharField(max_length=50)),
                ('alias', models.TextField()),
                ('source', models.TextField()),
            ],
            options={
                'db_table': 'chemical_aliases',
                'db_interaction_div': 'stitch',
            },
        ),
        migrations.CreateModel(
            name='ChemicalChemicalLinks',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('chemical1', models.CharField(max_length=50)),
                ('chemical2', models.CharField(max_length=50)),
                ('textmining', models.IntegerField()),
            ],
            options={
                'db_table': 'chemical_chemical_links',
                'db_interaction_div': 'stitch',
            },
        ),
        migrations.CreateModel(
            name='ChemicalChemicalLinksDetailed',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('chemical1', models.CharField(max_length=50)),
                ('chemical2', models.CharField(max_length=50)),
                ('similarity', models.IntegerField()),
                ('experimental', models.IntegerField()),
                ('database', models.IntegerField()),
                ('textmining', models.IntegerField()),
                ('combined_score', models.IntegerField()),
            ],
            options={
                'db_table': 'chemical_chemical_links_detailed',
                'db_interaction_div': 'stitch',
            },
        ),
        migrations.CreateModel(
            name='Chemicals',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('chemical_id', models.IntegerField()),
                ('preferred_name', models.TextField()),
                ('short_name', models.TextField()),
                ('mw', models.FloatField()),
                ('drug', models.BooleanField()),
            ],
            options={
                'db_table': 'chemicals',
                'db_interaction_div': 'stitch',
            },
        ),
        migrations.CreateModel(
            name='ChemicalsInchikeys',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('flat_chemical_id', models.CharField(max_length=50)),
                ('stereo_chemical_id', models.CharField(max_length=50)),
                ('source_cid', models.IntegerField()),
                ('inchikey', models.TextField()),
            ],
            options={
                'db_table': 'chemicals_inchikeys',
                'db_interaction_div': 'stitch',
            },
        ),
        migrations.CreateModel(
            name='ChemicalSources',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('flat_chemical', models.CharField(max_length=50)),
                ('stereo_chemical', models.CharField(max_length=50)),
                ('source_name', models.CharField(max_length=50)),
                ('source_id', models.CharField(max_length=50)),
            ],
            options={
                'db_table': 'chemical_sources',
                'db_interaction_div': 'stitch',
            },
        ),
        migrations.CreateModel(
            name='Collections',
            fields=[
                ('collection_id', models.IntegerField(serialize=False, primary_key=True)),
                ('pubmed_id', models.IntegerField()),
                ('comment', models.CharField(max_length=1000, blank=True)),
            ],
            options={
                'db_table': 'collections',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='EvidenceTransfers',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('target_protein_id_a', models.IntegerField()),
                ('target_protein_id_b', models.IntegerField()),
                ('source_protein_id_a', models.IntegerField()),
                ('source_protein_id_b', models.IntegerField()),
                ('transfer_score_c1', models.SmallIntegerField()),
                ('transfer_score_c2', models.SmallIntegerField()),
                ('source_ascore', models.SmallIntegerField()),
                ('source_escore', models.SmallIntegerField()),
                ('source_dscore', models.SmallIntegerField()),
                ('source_tscore', models.SmallIntegerField()),
            ],
            options={
                'db_table': 'evidence_transfers',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='Funccats',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('funccat_id', models.CharField(max_length=5)),
                ('funccat_description', models.CharField(max_length=100)),
            ],
            options={
                'db_table': 'funccats',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='FusionEvidence',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('target_protein_id_a', models.IntegerField()),
                ('target_protein_id_b', models.IntegerField()),
                ('source_protein', models.IntegerField()),
                ('source_species', models.IntegerField()),
                ('transfer_score_c1', models.SmallIntegerField()),
                ('transfer_score_c2', models.SmallIntegerField()),
                ('fusion_score', models.SmallIntegerField()),
            ],
            options={
                'db_table': 'fusion_evidence',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='Genes',
            fields=[
                ('gene_id', models.IntegerField(serialize=False, primary_key=True)),
                ('gene_external_id', models.CharField(max_length=100)),
                ('start_position_on_contig', models.IntegerField()),
                ('end_position_on_contig', models.IntegerField()),
                ('protein_size', models.IntegerField()),
            ],
            options={
                'db_table': 'genes',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='GenesProteins',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('protein_id', models.IntegerField()),
                ('gene_id', models.IntegerField()),
            ],
            options={
                'db_table': 'genes_proteins',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='ItemsAbstracts',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('protein_id', models.IntegerField()),
                ('abstract_id', models.CharField(max_length=20)),
                ('name', models.CharField(max_length=100)),
                ('abstract_length', models.IntegerField()),
                ('mesh_id', models.IntegerField()),
            ],
            options={
                'db_table': 'items_abstracts',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='Meshterms',
            fields=[
                ('mesh_id', models.IntegerField(serialize=False, primary_key=True)),
                ('description', models.CharField(max_length=255)),
            ],
            options={
                'db_table': 'meshterms',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='NodeNodeLinks',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('node_id_a', models.IntegerField()),
                ('node_type_b', models.IntegerField(null=True, blank=True)),
                ('node_id_b', models.IntegerField()),
                ('combined_score', models.SmallIntegerField(null=True, blank=True)),
                ('evidence_scores', fields.ArrayField(models.TextField())),
            ],
            options={
                'db_table': 'node_node_links',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='Orthgroups',
            fields=[
                ('orthgroup_id', models.IntegerField(serialize=False, primary_key=True)),
                ('orthgroup_external_id', models.CharField(max_length=20)),
                ('description', models.CharField(max_length=1000)),
                ('protein_count', models.IntegerField()),
                ('species_count', models.SmallIntegerField()),
            ],
            options={
                'db_table': 'orthgroups',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='OrthgroupsAbstracts',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('orthgroup_id', models.IntegerField()),
                ('abstract_id', models.CharField(max_length=250)),
            ],
            options={
                'db_table': 'orthgroups_abstracts',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='OrthgroupsFunccats',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('orthgroup_id', models.IntegerField()),
                ('funccat_id', models.IntegerField()),
            ],
            options={
                'db_table': 'orthgroups_funccats',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='OrthgroupsSets',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('orthgroup_id', models.IntegerField()),
                ('set_id', models.CharField(max_length=60)),
                ('is_database_set', models.BooleanField()),
            ],
            options={
                'db_table': 'orthgroups_sets',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='OrthgroupsSpecies',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('orthgroup_id', models.IntegerField()),
                ('species_id', models.IntegerField()),
                ('count', models.IntegerField()),
            ],
            options={
                'db_table': 'orthgroups_species',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='ProteinChemicalLinks',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('chemical', models.CharField(max_length=50)),
                ('protein', models.CharField(max_length=50)),
                ('combined_score', models.IntegerField()),
            ],
            options={
                'db_table': 'node_node_links',
                'db_interaction_div': 'stitch',
            },
        ),
        migrations.CreateModel(
            name='ProteinChemicalLinksDetailed',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('chemical', models.CharField(max_length=50)),
                ('protein', models.CharField(max_length=50)),
                ('experimental', models.IntegerField()),
                ('prediction', models.IntegerField()),
                ('database', models.IntegerField()),
                ('textmining', models.IntegerField()),
                ('combined_score', models.IntegerField()),
            ],
            options={
                'db_table': 'protein_chemical_links_detailed',
                'db_interaction_div': 'stitch',
            },
        ),
        migrations.CreateModel(
            name='ProteinChemicalLinksTransfer',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('chemical', models.CharField(max_length=50)),
                ('protein', models.CharField(max_length=50)),
                ('experimental_direct', models.IntegerField()),
                ('experimental_transferred', models.IntegerField()),
                ('prediction_direct', models.IntegerField()),
                ('prediction_transferred', models.IntegerField()),
                ('database_direct', models.IntegerField()),
                ('database_transferred', models.IntegerField()),
                ('textmining_direct', models.IntegerField()),
                ('textmining_transferred', models.IntegerField()),
                ('combined_score', models.IntegerField()),
            ],
            options={
                'db_table': 'protein_chemical_links_transfer',
                'db_interaction_div': 'stitch',
            },
        ),
        migrations.CreateModel(
            name='ProteinImageMatch',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('protein_id', models.IntegerField()),
                ('image_id', models.CharField(max_length=250)),
                ('identity', models.DecimalField(null=True, decimal_places=65535, max_digits=65535, blank=True)),
                ('source', models.CharField(max_length=10)),
                ('start_position_on_protein', models.IntegerField(null=True, blank=True)),
                ('end_position_on_protein', models.IntegerField(null=True, blank=True)),
                ('annotation', models.CharField(max_length=50, blank=True)),
            ],
            options={
                'db_table': 'protein_image_match',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='Proteins',
            fields=[
                ('protein_id', models.IntegerField(serialize=False, primary_key=True)),
                ('protein_external_id', models.CharField(max_length=100)),
                ('species_id', models.IntegerField()),
                ('protein_checksum', models.CharField(max_length=16)),
                ('protein_size', models.IntegerField()),
                ('annotation', models.CharField(max_length=600)),
                ('preferred_name', models.CharField(max_length=50)),
                ('annotation_word_vectors', models.TextField(blank=True)),
            ],
            options={
                'db_table': 'proteins',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='ProteinsMeshterms',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('mesh_id', models.IntegerField()),
                ('protein_id', models.IntegerField()),
            ],
            options={
                'db_table': 'proteins_meshterms',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='ProteinsNames',
            fields=[
                ('protein_name', models.CharField(max_length=60)),
                ('protein_id', models.IntegerField(serialize=False, primary_key=True)),
                ('species_id', models.IntegerField()),
                ('source', models.CharField(max_length=100)),
                ('is_preferred_name', models.NullBooleanField()),
                ('linkout', models.CharField(max_length=15, blank=True)),
            ],
            options={
                'db_table': 'proteins_names',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='ProteinsOrthgroups',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('orthgroup_id', models.IntegerField()),
                ('protein_id', models.IntegerField()),
                ('protein_external_id', models.CharField(max_length=50)),
                ('species_id', models.IntegerField()),
                ('start_position', models.IntegerField()),
                ('end_position', models.IntegerField()),
                ('preferred_name', models.CharField(max_length=50)),
                ('protein_annotation', models.CharField(max_length=100)),
                ('preferred_linkout_url', models.CharField(max_length=150, blank=True)),
            ],
            options={
                'db_table': 'proteins_orthgroups',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='ProteinsSequences',
            fields=[
                ('protein_id', models.IntegerField(serialize=False, primary_key=True)),
                ('sequence', models.TextField()),
            ],
            options={
                'db_table': 'proteins_sequences',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='ProteinsSmartlinkouts',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('protein_id', models.IntegerField()),
                ('protein_size', models.IntegerField()),
                ('smart_url', models.CharField(max_length=2000)),
            ],
            options={
                'db_table': 'proteins_smartlinkouts',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='Runs',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('run_id', models.IntegerField()),
                ('species_id', models.IntegerField()),
                ('contig_id', models.CharField(max_length=50)),
            ],
            options={
                'db_table': 'runs',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='RunsGenesProteins',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('run_id', models.IntegerField()),
                ('gene_id', models.IntegerField()),
                ('protein_id', models.IntegerField()),
                ('start_position_on_contig', models.IntegerField()),
                ('end_position_on_contig', models.IntegerField()),
                ('preferred_name', models.CharField(max_length=50)),
                ('annotation', models.CharField(max_length=100)),
            ],
            options={
                'db_table': 'runs_genes_proteins',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='RunsOrthgroups',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('run_id', models.IntegerField()),
                ('orthgroup_id', models.IntegerField()),
            ],
            options={
                'db_table': 'runs_orthgroups',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='ScoreTypes',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('score_id', models.SmallIntegerField()),
                ('score_type', models.CharField(max_length=35)),
            ],
            options={
                'db_table': 'score_types',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='Sets',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('set_id', models.CharField(max_length=60)),
                ('collection_id', models.IntegerField()),
                ('title', models.CharField(max_length=100, blank=True)),
                ('comment', models.CharField(max_length=255, blank=True)),
                ('url', models.CharField(max_length=255, blank=True)),
            ],
            options={
                'db_table': 'sets',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='SetsItems',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('set_id', models.IntegerField()),
                ('protein_id', models.IntegerField()),
                ('species_id', models.IntegerField()),
                ('is_database_set', models.BooleanField()),
            ],
            options={
                'db_table': 'sets_items',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='SetsPubmedrefs',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('set_id', models.CharField(max_length=60)),
                ('pubmed_id', models.CharField(max_length=20, blank=True)),
            ],
            options={
                'db_table': 'sets_pubmedrefs',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='Species',
            fields=[
                ('species_id', models.IntegerField(serialize=False, primary_key=True)),
                ('official_name', models.CharField(max_length=100)),
                ('compact_name', models.CharField(max_length=100)),
                ('kingdom', models.CharField(max_length=15)),
                ('type', models.CharField(max_length=10)),
            ],
            options={
                'db_table': 'species',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='SpeciesNames',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('species_id', models.IntegerField()),
                ('species_name', models.IntegerField()),
                ('official_name', models.CharField(max_length=255, blank=True)),
                ('is_string_species', models.NullBooleanField()),
            ],
            options={
                'db_table': 'species_names',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='SpeciesNodes',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('species_id', models.IntegerField()),
                ('species_name', models.CharField(max_length=255)),
                ('position', models.IntegerField()),
                ('size', models.IntegerField()),
            ],
            options={
                'db_table': 'species_nodes',
                'db_interaction_div': 'string',
            },
        ),
        migrations.CreateModel(
            name='Stitchactions',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('item_id_a', models.CharField(max_length=50)),
                ('item_id_b', models.CharField(max_length=50)),
                ('mode', models.CharField(max_length=50)),
                ('action', models.CharField(max_length=50, blank=True)),
                ('a_is_acting', models.IntegerField()),
                ('score', models.IntegerField()),
            ],
            options={
                'db_table': 'stitchactions',
                'db_interaction_div': 'stitch',
            },
        ),
        migrations.CreateModel(
            name='StitchNodeNode',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('node_id_a', models.IntegerField()),
                ('node_type_b', models.IntegerField()),
                ('node_id_b', models.IntegerField()),
                ('combined_score', models.SmallIntegerField()),
                ('evidence_scores', fields.ArrayField(models.TextField())),
            ],
            options={
                'db_table': 'node_node_links',
                'db_interaction_div': 'stitch',
            },
        ),
    ]
