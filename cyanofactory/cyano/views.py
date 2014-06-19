"""
Whole-cell knowledge base views

Author: Jonathan Karr, jkarr@stanford.edu
Affiliation: Covert Lab, Department of Bioengineering, Stanford University
Last updated: 2012-07-17

Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from __future__ import absolute_import

import os
from django.contrib.contenttypes.models import ContentType
from django.template.context import Context
from haystack.inputs import AutoQuery
import settings
import tempfile
from copy import deepcopy
from itertools import chain

from django.contrib.auth.decorators import login_required
from django.core.exceptions import ValidationError
from django.core.urlresolvers import reverse
from django.db import transaction
from django.db.models import Count, Sum
from django.db.models.fields import BooleanField, NullBooleanField, AutoField, BigIntegerField, DecimalField, FloatField, IntegerField, PositiveIntegerField, PositiveSmallIntegerField, SmallIntegerField
from django.db.models.fields.related import RelatedObject, ManyToManyField, ForeignKey
from django.db.models.query import EmptyQuerySet, QuerySet
from django.shortcuts import get_object_or_404
from django.utils.text import capfirst

from haystack.query import SearchQuerySet

from cyano.forms import ExportDataForm, ImportDataForm, ImportSpeciesForm
import cyano.helpers as chelpers
import cyano.models as cmodels
from cyano.models import PermissionEnum as perm
from cyano.decorators import resolve_to_objects, permission_required
from django.db.transaction import atomic
from django.http.response import HttpResponseRedirect, HttpResponseBadRequest, HttpResponse


def index(request):
    return chelpers.render_queryset_to_response(
        request,
        template = "cyano/index.html",
        )

@resolve_to_objects
@permission_required(perm.READ_NORMAL)
def species(request, species):
    contentCol1 = []
    contentCol2 = []
    contentCol3 = []

    if species is not None:
        contentCol1.append([
            [0, 'Compartments', cmodels.Compartment.objects.for_species(species).count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Compartment'})],
        ])
        
        chrs = cmodels.Chromosome.objects.for_species(species).values_list("length", "sequence")
        chrcontent = sum(chro[0] for chro in chrs)
        chr_gc_content = 0 if len(chrs) == 0 else sum([cmodels.Chromosome(sequence = chro[1]).get_gc_content() * chro[0] for chro in chrs]) / chrcontent

        plasmids = cmodels.Plasmid.objects.for_species(species).values_list("length", "sequence")
        plasmid_content = sum(plasmid[0] for plasmid in plasmids)
        plasmid_gc_content = 0 if len(plasmids) == 0 else sum([cmodels.Plasmid(sequence = plasmid[1]).get_gc_content() * plasmid[0] for plasmid in plasmids]) / plasmid_content

        contentCol1.append([
            [0, 'Genome', len(chrs) + len(plasmids), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Genome'})],
            [1, 'Chromosomes', len(chrs), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Chromosome'})],
            [2, 'Length', chrcontent, 'nt'],
            [2, 'GC-content', ('%0.1f' % (chr_gc_content * 100)), '%'],
            [1, 'Plasmids', len(plasmids), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Plasmid'})],
            [2, 'Length', plasmid_content, 'nt'],
            [2, 'GC-content', ('%0.1f' % (plasmid_gc_content * 100)), '%'],
        ])
                
        tus = cmodels.TranscriptionUnit.objects.for_species(species).annotate(num_genes = Count('genes'))
        nPolys = tus.filter(num_genes__gt = 1).count()
        contentCol1.append([
            [0, 'Transcription units', tus.count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'TranscriptionUnit'})],            
            [1, 'Monocistrons', tus.count() - nPolys],
            [1, 'Polycistrons', nPolys],
        ])
        
        genes = cmodels.Gene.objects.for_species(species).values_list("type__wid", "expression", "half_life")
        
        contentCol1.append([
            [0, 'Genes', len(genes), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Gene'})],
            [1, 'mRNA', sum((lambda x: x[0] == "mRNA")(x) for x in genes), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Gene'}) + '?type=mRNA'],
            [1, 'rRNA', sum((lambda x: x[0] == "rRNA")(x) for x in genes), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Gene'}) + '?type=rRNA'],
            [1, 'sRNA', sum((lambda x: x[0] == "sRNA")(x) for x in genes), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Gene'}) + '?type=sRNA'],
            [1, 'tRNA', sum((lambda x: x[0] == "tRNA")(x) for x in genes), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Gene'}) + '?type=tRNA'],
        ])
        
        chromosome_features = cmodels.ChromosomeFeature.objects.for_species(species).values_list("type__parent__wid", flat = True)
        chromosome_features_length = len(chromosome_features)
        chromosome_features_dnaa_box = sum((lambda x: x == "ChromosomeFeature-DnaA_box")(x) for x in chromosome_features)
        chromosome_features_str = sum((lambda x: x == "ChromosomeFeature-Short_Tandem_Repeat")(x) for x in chromosome_features)

        #contentCol1.append([
        #    [0, 'Chromosome features',
        #        chromosome_features_length,
        #        None,
        #        reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'ChromosomeFeature'}),
        #        ],
        #    [1, 'DnaA boxes',
        #        chromosome_features_dnaa_box,
        #        ],
        #    [1, 'Short tandem repeats',
        #        chromosome_features_str,
        #        ],
        #    [1, 'Other',
        #        chromosome_features_length - chromosome_features_dnaa_box - chromosome_features_str],
        #])
        
        metabolites = cmodels.Metabolite.objects.for_species(species).values_list("type__wid", "type__parent__wid", "type__parent__parent__wid", "biomass_composition", "media_composition")
        
        metabolites_aa_list = ["amino_acid", "modified_amino_acid", "non-standard_amino_acid", "vitamin_non-standard_amino_acid"]
        
        contentCol1.append([
            [0, 'Metabolites', len(metabolites), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Metabolite'})],
            [1, 'Amino acids', 
                sum((lambda x: x in metabolites_aa_list)(x) for x in metabolites)
                ],
            [1, 'Antibiotic',
                sum((lambda x: x[0] == "antibiotic")(x) for x in metabolites) +
                sum((lambda x: x[1] == "antibiotic")(x) for x in metabolites)
                ],
            [1, 'Gases', 
                sum((lambda x: x[0] == "gas")(x) for x in metabolites) +
                sum((lambda x: x[1] == "gas")(x) for x in metabolites)
                ],
            [1, 'Ions', 
                sum((lambda x: x[0] == "ion")(x) for x in metabolites) +
                sum((lambda x: x[1] == "ion")(x) for x in metabolites)
                ],
            [1, 'Lipids', 
                sum((lambda x: x[0] == "lipid")(x) for x in metabolites) +
                sum((lambda x: x[1] == "lipid")(x) for x in metabolites) +
                sum((lambda x: x[2] == "lipid")(x) for x in metabolites)
                ],
            [1, 'Vitamins',
                sum((lambda x: x[0] == "vitamin")(x) for x in metabolites) +
                sum((lambda x: x[1] == "vitamin")(x) for x in metabolites)
                ],
        ])

        mons = cmodels.ProteinMonomer.objects.for_species(species)
        mons_list = mons.values_list("localization__wid", "signal_sequence__type")
        mons_list_count = len(mons_list)
        monTermOrg_list = ["tc", "tm"]
        
        monDNABind = mons.filter(dna_footprint__length__gt=0).count()
        
        monLipo = sum((lambda x: x[1] == "lipoprotein")(x) for x in mons_list)
        monSecreted = sum((lambda x: x[1] == "secretory")(x) for x in mons_list)
        monTermOrg = sum((lambda x: x[0] in monTermOrg_list)(x) for x in mons_list)
        
        monIntMem = sum((lambda x: x[0] == "m" and x[1] != "lipoprotein")(x) for x in mons_list) +\
                    sum((lambda x: x[0] == "tm" and x[1] != "lipoprotein")(x) for x in mons_list)

        cpxs = cmodels.ProteinComplex.objects.for_species(species)
        cpxDNABind = cpxs.filter(dna_footprint__length__gt=0).count()
        cpxCount = cpxs.count()
        
        contentCol2.append([
            [0, 'Proteins', mons_list_count + cpxCount],
                [1, 'Monomers', mons_list_count, None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'ProteinMonomer'})],            
                    [2, 'DNA-binding', monDNABind],
                    [2, 'Integral membrane', monIntMem],
                    [2, 'Lipoprotein', monLipo, None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'ProteinMonomer'}) + '?signal_sequence__type=lipoprotein'],            
                    [2, 'Secreted', monSecreted, None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'ProteinMonomer'}) + '?signal_sequence__type=secretory'],
                    [2, 'Terminal organelle', monTermOrg],                    
                [1, 'Complexes', cpxCount, None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'ProteinComplex'})],
                    [2, 'DNA-binding', cpxDNABind],
        ])

        rxns = cmodels.Reaction.objects.for_species(species).values_list("processes__wid", flat = True)
        rxns_count = [
            len(rxns),
            sum((lambda x: x == "Process_DNADamage")(x) for x in rxns),
            sum((lambda x: x == "Process_DNARepair")(x) for x in rxns),
            sum((lambda x: x == "Process_Metabolism")(x) for x in rxns),
            sum((lambda x: x == "Process_ProteinDecay")(x) for x in rxns),
            sum((lambda x: x == "Process_ProteinModification")(x) for x in rxns),
            sum((lambda x: x == "Process_ReplicationInitiation")(x) for x in rxns),
            sum((lambda x: x == "Process_RNADecay")(x) for x in rxns),
            sum((lambda x: x == "Process_RNAModification")(x) for x in rxns),
            sum((lambda x: x == "Process_RNAProcessing")(x) for x in rxns),
            sum((lambda x: x == "Process_Transcription")(x) for x in rxns),
            sum((lambda x: x == "Process_Translation")(x) for x in rxns),
            sum((lambda x: x == "Process_tRNAAminoacylation")(x) for x in rxns)
        ]
        
        contentCol2.append([
            [0, 'Reactions', rxns_count[0], None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'})],
            [1, 'DNA damage', rxns_count[1], None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_DNADamage'],
            [1, 'DNA repair', rxns_count[2], None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_DNARepair'],
            [1, 'Metabolic', rxns_count[3], None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_Metabolism'],            
            [1, 'Protein decay', rxns_count[4], None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_ProteinDecay'],
            [1, 'Protein modification', rxns_count[5], None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_ProteinModification'],            
            [1, 'Replication Initiation', rxns_count[6], None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_ReplicationInitiation'],
            [1, 'RNA decay', rxns_count[7], None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_RNADecay'],
            [1, 'RNA modification', rxns_count[8], None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_RNAModification'],            
            [1, 'RNA processing', rxns_count[9], None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_RNAProcessing'],            
            [1, 'Transcription', rxns_count[10], None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_Transcription'],            
            [1, 'Translation', rxns_count[11], None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_Translation'],            
            [1, 'tRNA aminoacylation', rxns_count[12], None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_tRNAAminoacylation'],
            [1, 'Other', rxns_count[0] - sum(rxns_count[1:])],
        ])        
        
        tr = cmodels.TranscriptionalRegulation.objects.for_species(species).values_list("transcription_unit", "transcription_factor", "affinity", "activity")
        nTus = len(set([x[0] for x in tr]))
        nTfs = len(set([x[1] for x in tr]))
        contentCol3.append([
            [0, 'Transcriptional regulation'],
            [1, 'Interactions', tr.count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'TranscriptionalRegulation'})],
            [1, 'Transcriptional regulators', nTfs],
            [1, 'Regulated promoters', nTus],
        ])

        contentCol3.append([
            [0, 'Pathways', cmodels.Pathway.objects.filter(species__id = species.id).count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Pathway'})],
        ])
        
        stimuli = cmodels.Stimulus.objects.for_species(species).values_list("value", flat=True)
        nStimuli = sum((lambda x: x[3] is not None)(x) for x in stimuli)

        contentCol3.append([
            [0, 'Stimuli', len(stimuli), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Stimulus'})],
        ])
        
        nCellComp = sum((lambda x: x[3] is not None)(x) for x in metabolites)
        nMediaComp = sum((lambda x: x[4] is not None)(x) for x in metabolites)
        
        reactions = cmodels.Reaction.objects.for_species(species).\
            values_list("keq", "kinetics_forward__km", "kinetics_backward__km", "kinetics_forward__vmax", "kinetics_backward__vmax")
       
        nKineticsKeq = sum((lambda x: x[0] is not None)(x) for x in reactions)
        nKineticsKm = sum((lambda x: x[1] is not None)(x) for x in reactions) +\
                        sum((lambda x: x[2] is not None)(x) for x in reactions)
        nKineticsVmax = sum((lambda x: x[3] is not None)(x) for x in reactions) +\
                        sum((lambda x: x[4] is not None)(x) for x in reactions)

        nRnaExp = sum((lambda x: x[1] is not None)(x) for x in genes)
        nRnaHl = sum((lambda x: x[2] is not None)(x) for x in genes)

        nTrAffinity = sum((lambda x: x[2] is not None)(x) for x in tr)
        nTrActivity = sum((lambda x: x[3] is not None)(x) for x in tr)
       
        nOther = cmodels.Parameter.objects.for_species(species).count()
        nTotParameters = nCellComp + nMediaComp + nKineticsVmax + nRnaExp + nRnaHl + nStimuli + nTrAffinity + nTrActivity + nOther

        contentCol3.append([
            [0, 'Quantitative parameters', nTotParameters],
            [1, 'Cell composition', nCellComp],
            [1, 'Media composition', nMediaComp],            
            [1, 'Reaction K<sub>eq</sub>', nKineticsKeq],
            [1, 'Reaction K<sub>m</sub>', nKineticsKm],
            [1, 'Reaction V<sub>max</sub>', nKineticsVmax],
            [1, 'RNA expression', nRnaExp],
            [1, 'RNA half-lives', nRnaHl],
            [1, 'Stimulus values', nStimuli],   
            [1, 'Transcr. reg. activity', nTrAffinity],
            [1, 'Transcr. reg. affinity', nTrActivity],
            [1, 'Other', nOther, None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Parameter'})],            
        ])
        
        contentCol3.append([
            [0, 'Processes', cmodels.Process.objects.for_species(species).count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Process'})],
        ])
        
        contentCol3.append([
            [0, 'States', cmodels.State.objects.for_species(species).count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'State'})],
        ])

        contentCol3.append([
            [0, "Mass Spectrometry Data", cmodels.MassSpectrometryJob.objects.for_species(species).count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'MassSpectrometryJob'})],
        ])

    sources = {
        'total': 0,
        'types': [],
        'dates': [],
        'evidence_parameters': [],
        'evidence_species': [],
        'evidence_media': [],
        'evidence_pH': [],
        'evidence_temperature': [],
    }
    #if species is not None:
        #refs = cmodels.PublicationReference.objects.filter.for_species(species)
        #sources['total'] = refs.count()
        #sources['types'] = [
        #        {'type': 'Articles', 'count': refs.filter(type__wid='article').count()},
        #        {'type': 'Books', 'count': refs.filter(type__wid='book').count()},
        #        {'type': 'Thesis', 'count': refs.filter(type__wid='thesis').count()},
        #        {'type': 'Other', 'count': refs.filter(type__wid='misc').count()},
        #    ]
        #sources['dates'] = refs.filter(year__isnull=False).order_by('year').values('year').annotate(count=Count('year'))
        
            
        #nEstimated = cmodels.Parameter.objects.for_species(species).filter(value__evidence__is_experimentally_constrained=False).count()
        #nExpConstrained = nTotParameters - nEstimated
        #sources['evidence_parameters'] = [
        #     {'type': 'Experimentally constrained', 'count': nExpConstrained},
        #     {'type': 'Computationally estimated', 'count': nEstimated},
        #     ]
        
        # sources['evidence_species'] = cmodels.Evidence.objects.filter(species_component__species__id = species.id).values('species').annotate(count = Count('id'))
        # sources['evidence_media'] = cmodels.Evidence.objects.filter(species_component__species__id = species.id).values('media').annotate(count = Count('id'))
        # sources['evidence_pH'] = cmodels.Evidence.objects.filter(species_component__species__id = species.id).values('pH').annotate(count = Count('id'))
        # sources['evidence_temperature'] = cmodels.Evidence.objects.filter(species_component__species__id = species.id).values('temperature').annotate(count = Count('id'))

    printCol1 = []
    printCol2 = []
    printCol3 = []

    i = 0
    for x in contentCol1:
        i += 1
        for y in x:
            printCol1.append([i] + y)
    i = 0
    for x in contentCol2:
        i += 1
        for y in x:
            printCol2.append([i] + y)
    i = 0
    for x in contentCol3:
        i += 1
        for y in x:
            printCol3.append([i] + y)

    return chelpers.render_queryset_to_response(
        species = species,
        data = {
            'content': [printCol1, printCol2, printCol3],
            'contentRows': range(max(len(printCol1), len(printCol2), len(printCol3))),
            'sources': sources,    
            },
        request = request, 
        template = 'cyano/species.html')

@resolve_to_objects
def about(request, species=None):
    return chelpers.render_queryset_to_response(
        species = species,
        request = request, 
        template = 'cyano/about.html', 
        data = {
            'ROOT_URL': settings.ROOT_URL,
        }
    )
        
@resolve_to_objects
def tutorial(request, species=None):
    return chelpers.render_queryset_to_response(
        species = species,
        request = request, 
        template = 'cyano/tutorial.html')

@resolve_to_objects
def licensing(request, species=None):
    return chelpers.render_queryset_to_response(
        species = species,
        request = request, 
        template = 'cyano/license.html')

@login_required
@resolve_to_objects
def users(request, species = None):
    queryset = cmodels.UserProfile.objects.all().filter(user__is_active = True)
    return chelpers.render_queryset_to_response(
        species = species,
        request = request, 
        models = [cmodels.UserProfile],
        queryset = queryset,
        template = 'cyano/users.html')
        
@login_required   
@resolve_to_objects
def user(request, username, species = None):
    queryset = chelpers.objectToQuerySet(get_object_or_404(cmodels.UserProfile, user__username = username), model = cmodels.UserProfile)
    return chelpers.render_queryset_to_response(
        species = species,
        request = request,
        models = [cmodels.UserProfile],
        queryset = queryset,
        template = 'cyano/user.html')    

@resolve_to_objects
def search(request, species = None):
    query = request.GET.get('q', '')
    engine = request.GET.get('engine', 'haystack')
    
    if engine == 'haystack' or not getattr(settings, 'GOOGLE_SEARCH_ENABLED', False):
        return search_haystack(request, species, query)
    else:
        return search_google(request, species, query)
                
def search_haystack(request, species, query):
    results = SearchQuerySet().filter(species_wid=species.wid).filter(content=AutoQuery(query))

    #calculate facets      
    facets = results.facet('model_type')
    ##print facets
    ##print facets.facet_counts()

    models = []
    model_name_facet = []

    if request.is_ajax():
        template = "cyano/search_page.html"
    else:
        template = "cyano/search.html"

    model_type = request.GET.get('model_type', '')

    if facets.facet_counts():
        tmp = facets.facet_counts()['fields']['model_type']
        ##print tmp
        for tmp2 in tmp:
            ##print "tmp2", tmp2
            model_name = cmodels.TableMeta.objects.get(model_name__iexact=tmp2[0]).model_name
            model_name_facet.append({
                'name':model_name,
                'verbose_name': chelpers.getModel(model_name)._meta.verbose_name,
                'count':tmp2[1],
                })
            models.append(chelpers.getModel(model_name))
        model_name_facet.sort(lambda x, y:cmp(x['verbose_name'], y['verbose_name']))

        #narrow search by facets
        if model_type:
            results = results.filter(model_type=model_type)

    #order results
    results = results.order_by('wid')

    results.model = cmodels.Entry

    #for result in results:
    #    print result.model_name, result
    
    #form response
    return chelpers.render_queryset_to_response(
        species = species,
        request = request, 
        models = models,
        queryset = results,
        template = template,
        data = {
            'query': query,
            'engine': 'haystack',
            'model_type': model_type,
            'modelNameFacet': model_name_facet,
            })

def search_google(request, species, query):
    return chelpers.render_queryset_to_response(
        species = species,
        request = request, 
        template = 'cyano/googleSearch.html', 
        data = {
            'query': query,
            'engine': 'google',
            })

@resolve_to_objects
@permission_required(perm.READ_NORMAL)
def listing(request, species, model):
    from itertools import groupby
    from collections import OrderedDict

    objects = model.objects.for_species(species)

    facet_fields = []    
    for field_full_name in model._meta.facet_fields:
        #facet
        field_names = str(field_full_name).split('__')
        tmp_model = model
        field_verbose_name = []
        for field_name in field_names:
            field = tmp_model._meta.get_field_by_name(field_name)[0]
            field_verbose_name.append(field.verbose_name)
            if isinstance(field, (ForeignKey, ManyToManyField)):
                tmp_model = field.rel.to
        field_verbose_name = ' &#8250; '.join(field_verbose_name)
                
        if isinstance(field, (ForeignKey, ManyToManyField)) and not issubclass(field.rel.to, cmodels.Entry):
            continue
        
        if isinstance(field, (ForeignKey, ManyToManyField)):
            tmp = model.objects.for_species(species).order_by(field_full_name + '__name').values(field_full_name).annotate(count=Count(field_full_name))
        else:
            tmp = model.objects.for_species(species).order_by(field_full_name).values(field_full_name).annotate(count=Count(field_full_name))
        facets = []
        for facet in tmp:
            value = facet[field_full_name]            
            if value is None or unicode(value) == '':
                continue
            
            if isinstance(field, (ForeignKey, ManyToManyField)):
                tmp2 = tmp_model.objects.values('wid', 'name').get(id=value)
                id_ = tmp2['wid']
                name = capfirst(tmp2['name'])
            elif (field.choices is not None) and (len(field.choices) > 0) and (not isinstance(field, (BooleanField, NullBooleanField))):    
                id_ = value
                choices = [choice[0] for choice in field.choices]
                if id_ in choices:
                    name = field.choices[choices.index(id_)][1]
                else:
                    name = capfirst(value)
            else:
                id_ = value
                name = capfirst(value)
            if value is not None and unicode(value) != '':
                facets.append({
                    'id': unicode(id_), 
                    'name': unicode(name or id_),
                    'count': facet['count']})
        if len(facets) > 1:
            facet_fields.append({ 
                'name': field_full_name,
                'verbose_name': field_verbose_name, 
                'facets': facets,
                })
    
        #filter
        val = request.GET.get(field_full_name)
        if val:
            if isinstance(field, (ForeignKey, ManyToManyField)):
                kwargs = {field_full_name + '__wid': val}
            elif isinstance(field, (BooleanField, NullBooleanField)):
                kwargs = {field_full_name: val == 'True'}
            elif isinstance(field, (AutoField, BigIntegerField, DecimalField, FloatField, IntegerField, PositiveIntegerField, PositiveSmallIntegerField, SmallIntegerField)):
                kwargs = {field_full_name: float(val)}
            else:
                kwargs = {field_full_name: val}
            objects = objects.filter(**kwargs)

    if request.is_ajax():
        template = "cyano/list_page.html"
    else:
        template = "cyano/list.html"

    groups = None

    if hasattr(model._meta, "group_field"):
        field_name = getattr(model._meta, "group_field")
        group_field = model._meta.get_field_by_name(getattr(model._meta, "group_field"))[0]

        if group_field:
            objects = objects.order_by(field_name, "wid")

            #if isinstance(group_field, (ForeignKey, ManyToManyField)):
            #    print "m2m or fk"
            #else:
            #    print "other"

            objects = objects.prefetch_related(group_field.name)

            def group_func(x):
                try:
                    return getattr(x, field_name).all()[0].wid
                except IndexError:
                    return "None"

            groups = OrderedDict((k, list(v)) for k, v in groupby(objects, group_func))
            #print groups
    else:
        objects = objects.order_by("wid")
    
    return chelpers.render_queryset_to_response(
        species=species,
        request=request,
        models=[model],
        queryset=objects,
        template=template,
        data={
            'groups': groups,
            'facet_fields': facet_fields
        }
    )

@resolve_to_objects
@permission_required(perm.READ_NORMAL)
def detail(request, species, model, item):
    fieldsets = deepcopy(model._meta.fieldsets)
    
    if request.GET.get('format', 'html') == "html":
        #filter out type, metadata
        fieldset_names = [x[0] for x in filter(lambda x: isinstance(x, tuple), fieldsets)]
        if 'Type' in fieldset_names:
            idx = fieldset_names.index('Type')
            del fieldsets[idx]
            
        #filter out empty fields
        fieldsets = chelpers.create_detail_fieldset(species, item, fieldsets, request.user.is_anonymous())

    qs = chelpers.objectToQuerySet(item, model = model)

    #render response
    return chelpers.render_queryset_to_response(
        species = species,        
        request = request, 
        models = [model],
        queryset = qs,
        template = 'cyano/detail.html', 
        data = {
            'fieldsets': fieldsets,
            'message': request.GET.get('message', ''),
            })


@resolve_to_objects
@permission_required(perm.READ_NORMAL)
def detail_field(request, species, model, item):
    from django.template import loader
    from django.utils.html import strip_tags

    if request.GET.get('name') is None:
        return HttpResponseBadRequest("Unknown field")

    strip = request.GET.get('strip', False)

    output = chelpers.format_field_detail_view(species, item, request.GET.get('name'), request.user.is_anonymous())

    if output is None:
        return HttpResponseBadRequest("Unknown field")

    template = loader.get_template("cyano/field.html")
    c = Context({'request': request, 'data': output})

    rendered = template.render(c)

    if strip:
        rendered = strip_tags(rendered)

    return HttpResponse(rendered)

@resolve_to_objects
@permission_required(perm.READ_HISTORY)
def history(request, species, model=None, item=None):
    revisions = []
    entry = []
    date = None

    if item:
        # Item specific
        objects = cmodels.Revision.objects.filter(object_id=item.pk).distinct().order_by("-detail")

    elif model:
        # Model specific
        components = model.objects.for_species(species)
        ct_id = ContentType.objects.get_for_model(model).pk
        objects = cmodels.Revision.objects.filter(object_id__in=components, content_type__pk=ct_id).order_by("-detail")

    else:
        # Whole species specific
        components = cmodels.SpeciesComponent.objects.for_species(species)
        objects = cmodels.Revision.objects.filter(object_id__in=components).distinct().order_by("-detail")
    
    for obj in objects:
        if not issubclass(ContentType.objects.get_for_id(obj.content_type_id).model_class(), cmodels.Entry):
            continue

        last_date = date
        wid = obj.current.wid
        item_model = obj.current.model_type.model_name
        detail_id = obj.detail.pk
        date = obj.detail.date.date()
        time = obj.detail.date.strftime("%H:%M")
        reason = obj.detail.reason
        author = obj.detail.user
        url = reverse("cyano.views.history_detail", kwargs = {"species_wid": species.wid, "model_type": item_model, "wid": wid, "detail_id": detail_id})

        if last_date != date:
            revisions.append(entry)
            entry = [date, []]

        entry[1].append({'id': detail_id, 'time': time, 'wid': wid, 'reason': reason, 'author': author, 'url': url})
    revisions.append(entry)

    if item:
        qs = chelpers.objectToQuerySet(item, model = model)
    else:
        qs = objects

    return chelpers.render_queryset_to_response(
        species=species,
        request=request,
        models=[model],
        queryset=qs,
        template='cyano/history.html',
        data = {
            'revisions': revisions
            })

@resolve_to_objects
@permission_required(perm.READ_HISTORY)
def history_detail(request, species, model, item, detail_id):
    fieldsets = deepcopy(model._meta.fieldsets)
    
    #filter out type, metadata
    fieldset_names = [x[0] for x in fieldsets]
    if 'Type' in fieldset_names:
        idx = fieldset_names.index('Type')
        del fieldsets[idx]
        
    #filter out empty fields
    item = chelpers.get_history(species, item, detail_id)

    rmfieldsets = []
    for idx in range(len(fieldsets)):
        rmfields = []
        for idx2 in range(len(fieldsets[idx][1]['fields'])):
            if isinstance(fieldsets[idx][1]['fields'][idx2], dict):
                field_name = fieldsets[idx][1]['fields'][idx2]['name']
                verbose_name = fieldsets[idx][1]['fields'][idx2]['verbose_name']
            else:
                field_name = fieldsets[idx][1]['fields'][idx2]
                field = model._meta.get_field_by_name(field_name)[0]
                if isinstance(field, RelatedObject):
                    verbose_name = capfirst(field.get_accessor_name())
                else:
                    verbose_name = field.verbose_name
                
            data = chelpers.format_field_detail_view(species, item, field_name, request.user.is_anonymous(), detail_id)
            if (data is None) or (data == ''):
                rmfields = [idx2] + rmfields
            
            fieldsets[idx][1]['fields'][idx2] = {'verbose_name': verbose_name.replace(" ", '&nbsp;').replace("-", "&#8209;"), 'data': data}
        for idx2 in rmfields:
            del fieldsets[idx][1]['fields'][idx2]
        if len(fieldsets[idx][1]['fields']) == 0:
            rmfieldsets = [idx] + rmfieldsets
    for idx in rmfieldsets:
        del fieldsets[idx]
    
    #form query set
    qs = chelpers.objectToQuerySet(item, model=model)
    ct_id = ContentType.objects.get_for_model(model).pk

    # prev rev
    prev_rev = cmodels.Revision.objects.filter(object_id=item.pk, content_type__pk=ct_id, detail_id__lt=detail_id).distinct().order_by("-detail").first()
    if prev_rev:
        prev_rev = prev_rev.detail

    new_rev = cmodels.Revision.objects.filter(object_id=item.pk, content_type__pk=ct_id, detail_id__gt=detail_id).distinct().order_by("detail").first()
    if new_rev:
        new_rev = new_rev.detail

    latest_rev = cmodels.Revision.objects.filter(object_id=item.pk, content_type__pk=ct_id).distinct().order_by("-detail").first().detail

    #render response
    return chelpers.render_queryset_to_response(
        species=species,
        request=request,
        models=[model],
        queryset=qs,
        template='cyano/history_detail.html',
        data={
            'fieldsets': fieldsets,
            'message': request.GET.get('message', ''),
            'latest_revision': latest_rev,
            'previous_revision': prev_rev,
            'revision': cmodels.RevisionDetail.objects.get(pk=detail_id),
            'newer_revision': new_rev
        }
    )

@login_required
@resolve_to_objects
@permission_required(perm.WRITE_NORMAL)
def add(request, species=None, model=None):
    return edit(request, species=species, model=model, action='add')

@login_required
@resolve_to_objects 
@permission_required(perm.WRITE_NORMAL)
def edit(request, species, model = None, item = None, action='edit'):
    from collections import defaultdict
    
    #retrieve object
    if action == 'add':
        obj = model()
    else:
        if item is None:
            obj = species
        else:
            obj = item
        
        if model is None:
            model = obj.__class__
    
    #save object
    error_messages = {}
    if request.method == 'POST':
            submitted_data = chelpers.get_edit_form_data(model, request.POST, user=request.user.profile)
            
            data = submitted_data
            data['id'] = obj.id
            data['species'] = species.wid
            data['model_type'] = model.__name__
            data['wid'] = data['wid'] if action == 'add' else obj.wid
            
            try:
                with transaction.atomic():
                    #validate is WID unique
                    if issubclass(model, cmodels.SpeciesComponent):
                        qs = cmodels.SpeciesComponent.objects.values('wid', 'model_type__model_name').filter(species__wid=species.wid)
                    else:
                        qs = model.objects.values('wid', 'model_type__model_name').all()

                    if action == 'edit':
                        qs = qs.exclude(id=obj.id)

                    wids = defaultdict(list)
                    for x in qs:
                        wids[x['wid']].append(x['model_type__model_name'])

                    if data['wid'] in wids.keys():
                        if model.__name__ in wids[data['wid']]:
                            raise ValidationError({'wid': 'Value must be unique for model'})

                    wids[data['wid']] = model.__name__

                    #validate
                    data = chelpers.validate_object_fields(model, data, wids, species.wid, data['wid'])
                    chelpers.validate_revision_detail(data)
                    chelpers.validate_model_objects(model, data)
                    chelpers.validate_model_unique(model, [data])

                    #save
                    obj = chelpers.save_object_data(species, obj, data, {}, request.user, save=False, save_m2m=False)
                    obj = chelpers.save_object_data(species, obj, data, {data['wid']: obj}, request.user, save=True, save_m2m=False)
                    obj = chelpers.save_object_data(species, obj, data, {data['wid']: obj}, request.user, save=True, save_m2m=True)

                    #redirect to details page
                    return HttpResponseRedirect(obj.get_absolute_url(species))
            except ValidationError as error:
                if hasattr(error, "message_dict"):
                    error_messages = error.message_dict
                else:
                    raise

    #form query set
    if action == 'edit':
        qs = chelpers.objectToQuerySet(obj, model = model)
    else:
        obj = None
        qs = model.objects.none()
        
    #display form
    fields, initial_values = chelpers.get_edit_form_fields(None if species is None else species.wid, model, obj=obj)
    
    if request.method == 'POST':
        initial_values = submitted_data
    return chelpers.render_queryset_to_response(
        species = species,
        request = request, 
        models = [model],
        queryset = qs,
        template = 'cyano/edit.html', 
        data = {
            'action': action,
            'fields': fields,
            #'references_choices': cmodels.PublicationReference.objects.filter(species__wid = species.wid).values_list('wid'),
            'initial_values': initial_values,
            'error_messages': error_messages,
            }
        )

@login_required
@resolve_to_objects 
@permission_required(perm.WRITE_DELETE)
def delete(request, species, model = None, item = None):
    #retrieve object
    if item is None:
        obj = species
    else:
        obj = item
        
    if model is None:
        model = obj.__class__
    
    qs = chelpers.objectToQuerySet(obj, model = model)
    
    #delete
    if request.method == 'POST':
        # Todo: Should be revisioned with custom message
        rev_detail = cmodels.RevisionDetail(user=request.user.profile, reason="Delete "+item.wid)
        obj.delete(species, rev_detail)
        return HttpResponseRedirect(reverse('cyano.views.listing', kwargs={'species_wid':species.wid, 'model_type': model.__name__}))
        
    #confirmation message
    return chelpers.render_queryset_to_response(
        species = species,
        request = request, 
        models = [model],
        queryset = qs,
        template = 'cyano/delete.html', 
        )

@resolve_to_objects
@permission_required(perm.READ_NORMAL)
def exportData(request, species):
    form = ExportDataForm(None)
    if not form.is_valid():        
        return chelpers.render_queryset_to_response(
            species=species,
            request = request,
            template = 'cyano/exportDataForm.html', 
            data = {
                'form': form
                }
            )
    else:
        species = species.objects.get(wid = form.cleaned_data['species'])
        queryset = EmptyQuerySet()
        models = []
        
        if form.cleaned_data['all_model_types'] == 'True':
            model_types = chelpers.getObjectTypes()
        else:
            model_types = form.cleaned_data['model_type']
        
        for model_type in model_types:
            model = chelpers.getModel(model_type)
            if issubclass(model, cmodels.SpeciesComponent):
                queryset = chain(queryset, model.objects.for_species(species).select_related(depth=2).all())
            else:
                queryset = chain(queryset, model.objects.for_species(species).select_related(depth=2))
            models.append(chelpers.getModel(model_type))
        
        return chelpers.render_queryset_to_response(
            species = species,
            request = request, 
            queryset = queryset, 
            template = 'cyano/exportDataResult.html', 
            models = models)

@login_required
@resolve_to_objects
@permission_required(perm.WRITE_NORMAL)
def importData(request, species=None):
    data = {}
    
    if request.method == 'POST':
        form = ImportDataForm(request.POST, request.FILES)
        
        if form.is_valid():
            if form.cleaned_data.get('species'):
                selected_species_wid = form.cleaned_data.get('species')
            else:
                selected_species_wid = species.wid
            
            if selected_species_wid:
                #save to temporary file
                #originalFileName, originalFileExtension = os.path.splitext(request.FILES['file'].name)[1]
                originalFileExtension = os.path.splitext(request.FILES['file'].name)[1]
                fid = tempfile.NamedTemporaryFile(suffix = originalFileExtension, delete = False)
                filename = fid.name
                for chunk in request.FILES['file'].chunks():
                    fid.write(chunk)
                fid.close()
                
                #read file
                data_type = form.cleaned_data["data_type"]
                
                args = {"filename": filename,
                        "wid": selected_species_wid,
                        "user": request.user.username,
                        "reason": form.cleaned_data['reason']}

                if data_type == "genbank":
                    from cyano.tasks import genbank
                    genbank.delay(name = form.cleaned_data["chromosome"],
                                  chromosome = form.cleaned_data["chromosome_wid"],
                                  **args)
                elif data_type == "sbml":
                    from cyano.tasks import sbml
                    sbml.delay(**args)
                elif data_type == "proopdb":
                    from cyano.tasks import proopdb
                    proopdb.delay(**args)
                
                data['success'] = 'success'
                data['message'] = "New import job created for %s" % (request.FILES['file'].name)

    else:
        form = ImportDataForm(None)

    data["form"] = form
    return chelpers.render_queryset_to_response(
        species = species,
        request = request, 
        template = 'cyano/importDataForm.html', 
        data = data
        )

@login_required
@resolve_to_objects
@permission_required(perm.WRITE_NORMAL)
@atomic
def importSpeciesData(request, species=None):
    data = {}
    
    if request.method == 'POST':
        form = ImportSpeciesForm(request.POST)
        
        if form.is_valid():            
            rev = cmodels.RevisionDetail(user = request.user.profile, reason = form.cleaned_data["reason"])
            rev.save()
            
            mutant = False
            if species: # Create mutant
                import itertools
                
                mutant = True
                
                through = cmodels.SpeciesComponent.species.through
                    
                component = cmodels.SpeciesComponent.objects.for_species(species).values_list("pk", flat=True).order_by("pk")
                
                species.pk = None
                species.id = None
                species.wid = form.cleaned_data['new_wid']
                species.name = form.cleaned_data['new_species']
                species.save(rev)
                
                through.objects.bulk_create(map(lambda x: through(species_id = x[0], speciescomponent_id = x[1]), itertools.izip(itertools.cycle([species.pk]), component)))
            else: # Create species
                species = cmodels.Species(wid = form.cleaned_data['new_wid'], name = form.cleaned_data['new_species'])
                species.save(rev)
                cmodels.Pathway.add_boehringer_pathway(species, rev)

            # Assign permissions
            new_perm, _ = cmodels.UserPermission.objects.get_or_create(entry = species, user = request.user.profile)
            new_perm.allow.add(*cmodels.Permission.objects.all())

            data['success'] = 'success'
            data['message'] = "New %s %s created" % ("mutant" if mutant else "species", species.name)

    else:
        form = ImportSpeciesForm(None)

    data["form"] = form
    return chelpers.render_queryset_to_response(
        species = species,
        request = request, 
        template = 'cyano/importSpeciesForm.html', 
        data = data,
        )
    
def validate(request, species_wid):
    errors = chelpers.get_invalid_objects(cmodels.Species.objects.values('id').get(wid=species_wid)['id'])
    
    return chelpers.render_queryset_to_response(
        species = species,
        request = request, 
        template = 'cyano/validate.html', 
        data = {            
            'errors': errors
            },
        )

@login_required
@resolve_to_objects
def password_change_required(request, species=None):
    from django.http import HttpResponseRedirect
    from django.contrib.auth.views import password_change
    from django.contrib.auth.forms import AdminPasswordChangeForm
    
    if not request.user.profile.force_password_change:
        return HttpResponseRedirect(reverse("cyano.views.index"))
    
    context = chelpers.get_extra_context(
        species = species,
        request = request)
    
    return password_change(request,
                           "registration/password_change_form_required.html",
                           password_change_form = AdminPasswordChangeForm,
                           extra_context = context)

@resolve_to_objects
def login(request, species=None):
    from django.contrib.auth.views import login as djlogin
    from urllib import unquote

    msg = request.GET.get("message", "")

    context = chelpers.get_extra_context(
        species=species,
        request=request,
    )

    context['message'] = unquote(msg)[:50]

    return djlogin(request, extra_context = context)

@resolve_to_objects
def logout(request, species=None):
    from django.contrib.auth.views import logout as djlogout

    context = chelpers.get_extra_context(
        species = species,
        request = request)
    
    return djlogout(request, extra_context = context)

def sitemap(request):
    return chelpers.render_queryset_to_response(
        request = request, 
        template = 'cyano/sitemap.xml', 
        data = {
            'ROOT_URL': settings.ROOT_URL,
            'qs_species': cmodels.Species.objects.all(),
        }
    )
    
def sitemap_toplevel(request):
    return chelpers.render_queryset_to_response(
        request = request, 
        template = 'cyano/sitemap_toplevel.xml', 
        data = {
            'ROOT_URL': settings.ROOT_URL,
        }
    )

@resolve_to_objects
def sitemap_species(request, species):
    return chelpers.render_queryset_to_response(
        species = species,
        request = request, 
        template = 'cyano/sitemap_species.xml',
        data = {
            'ROOT_URL': settings.ROOT_URL,
            'entries': cmodels.SpeciesComponent.objects.for_species(species).select_related("detail"),
        }
    )

@login_required
@resolve_to_objects
@permission_required(perm.WRITE_PERMISSION)
def permission_edit(request, species, model = None, item = None):
    return permission(request, species = species, model = model, item = item, edit = True)

@resolve_to_objects
@permission_required(perm.READ_PERMISSION)
def permission(request, species, model = None, item = None, edit = False):
    users = cmodels.UserProfile.objects.all().filter(user__is_active = True)
    groups = cmodels.GroupProfile.objects.all()
    
    # Permissions for item or species depending on page
    entry = species if item == None else item

    if edit:
        if request.method == 'POST':
            # Form submit -> Save changes
            post_types = ["uid_allow", "gid_allow", "uid_deny", "gid_deny"]
            error = False
            
            # Check that all needed fields are in the request
            if not all(p in request.POST.keys() for p in post_types):
                error = True
    
            if not error:
                new_permissions = {}
                
                # Make QueryDict writable
                request.POST = request.POST.copy()
                
                # Add missing fields to the POST request
                for k in cmodels.Permission.permission_types + post_types:
                    if not k in request.POST:
                        request.POST[k] = None
                
                perm_len = None
    
                # Remove all fields except the needed ones, verify
                # that all new permissions are numbers
                for k, v in request.POST.lists():                                
                    if k in cmodels.Permission.permission_types + post_types:
                        try:
                            # was missing
                            if v[0] == None:
                                new_permissions[k] = []
                            else:
                                new_permissions[k] = map(lambda x: int(x), v)
    
                            if k in cmodels.Permission.permission_types:
                                # Verify that all permission fields have the same length
                                if perm_len == None:
                                    perm_len = len(new_permissions[k])
                                else:
                                    
                                    if perm_len != len(new_permissions[k]):
                                        raise ValueError()
                        except ValueError:
                            error = True
    
            if not error:
                # Remove id 0 caused by "Add new" field
                new_permissions["uid_allow"] = new_permissions["uid_allow"][:-1]
                new_permissions["gid_allow"] = new_permissions["gid_allow"][:-1]
                new_permissions["uid_deny"] = new_permissions["uid_deny"][:-1]
                new_permissions["gid_deny"] = new_permissions["gid_deny"][:-1]
                
                user_count_allow = len(new_permissions["uid_allow"])
                group_count_allow = len(new_permissions["gid_allow"])
                user_count_deny = len(new_permissions["uid_deny"])
                group_count_deny = len(new_permissions["gid_deny"])
    
                # Verify that user and group ids are valid
                user_verify = cmodels.UserProfile.objects.filter(pk__in = new_permissions["uid_allow"]).count()
                user_verify += cmodels.UserProfile.objects.filter(pk__in = new_permissions["uid_deny"]).count()
                group_verify = cmodels.GroupProfile.objects.filter(pk__in = new_permissions["gid_allow"]).count()
                group_verify += cmodels.GroupProfile.objects.filter(pk__in = new_permissions["gid_deny"]).count()
                
                # If number of permissions and users+groups mismatches
                new_perm_count = user_count_allow + group_count_allow + user_count_deny + group_count_deny
                
                if len(new_permissions["FULL_ACCESS"]) != new_perm_count:
                    error = True
    
                # If number of submitted data larger then database data at least one item was invalid
                if new_perm_count != user_verify + group_verify:
                    error = True
    
            if error:
                return chelpers.render_queryset_to_response_error(
                    request,
                    species = species,
                    error = 400,
                    msg = "Invalid request on permission page.",
                    msg_debug = "POST was: " + str(request.POST.lists()))        
            
            # Data is valid, begin with database operations
            
            new_permissions_user_allow = {}
            new_permissions_group_allow = {}
            new_permissions_user_deny = {}
            new_permissions_group_deny = {}
    
            # Rotate array to user -> permissions
            for i, user in enumerate(new_permissions["uid_allow"]):
                new_permissions_user_allow[cmodels.UserProfile.objects.get(pk = user)] = [new_permissions[x][i] for x in cmodels.Permission.permission_types]
            
            # Rotate array to group -> permissions
            for i, group in enumerate(new_permissions["gid_allow"], user_count_allow):  
                new_permissions_group_allow[cmodels.GroupProfile.objects.get(pk = group)] = [new_permissions[x][i] for x in cmodels.Permission.permission_types]
            
            # Rotate array to user -> permissions
            for i, user in enumerate(new_permissions["uid_deny"], user_count_allow + group_count_allow):
                new_permissions_user_deny[cmodels.UserProfile.objects.get(pk = user)] = [new_permissions[x][i] for x in cmodels.Permission.permission_types]
            
            # Rotate array to group -> permissions
            for i, group in enumerate(new_permissions["gid_deny"], user_count_allow + group_count_allow + user_count_deny):  
                new_permissions_group_deny[cmodels.GroupProfile.objects.get(pk = group)] = [new_permissions[x][i] for x in cmodels.Permission.permission_types]

            # Test if user is missing but had permission before (Happens with "Remove this user")
            all_user_perms = cmodels.UserPermission.objects.filter(entry = entry).prefetch_related("user")
            for permission in all_user_perms:
                if not permission.user in new_permissions_user_allow and not permission.user in new_permissions_user_deny:
                    permission.delete()
                else:
                    if not permission.user in new_permissions_user_allow:
                        permission.allow.clear()
                    if not permission.user in new_permissions_user_deny:
                        permission.deny.clear()
            
            # Test if group is missing but had permission before (Happens with "Remove this group")
            all_group_perms = cmodels.GroupPermission.objects.filter(entry = entry).prefetch_related("group")
            for permission in all_group_perms:
                if not permission.group in new_permissions_group_allow and not permission.group in new_permissions_group_deny:
                    permission.delete()
                else:
                    if not permission.group in new_permissions_group_allow:
                        permission.allow.clear()
                    if not permission.group in new_permissions_group_deny:
                        permission.deny.clear()

            # Create or alter the permission object associated to entry+user/group depending on the
            # new permission settings
            for u, p in new_permissions_user_allow.iteritems():
                u._handle_permission_list(entry, p, allow=True)

            for u, p in new_permissions_user_deny.iteritems():
                u._handle_permission_list(entry, p, allow=False)
                        
            for g, p in new_permissions_group_allow.iteritems():
                g._handle_permission_list(entry, p, allow=True)
                
            for g, p in new_permissions_group_deny.iteritems():
                g._handle_permission_list(entry, p, allow=False)

    # Rendering of the permission page
    user_permissions_db = cmodels.UserPermission.objects.filter(entry = entry).prefetch_related("allow", "deny", "user", "user__user")
    group_permissions_db = cmodels.GroupPermission.objects.filter(entry = entry).prefetch_related("allow", "deny", "group", "group__group")

    user_permissions_allow = []
    group_permissions_allow = []
    user_permissions_deny = []
    group_permissions_deny = []
    
    # Calculate permission array for user and groups
    # Contains user mapping to 8 numbers, 1 if permission available, 0 if not
    for permission in user_permissions_db:
        user_perm_allow = [permission.user.id, [1 if cmodels.Permission.get_by_pk(x) in permission.allow.all() else 0 for x in range(1, 9)]]
        user_perm_deny = [permission.user.id, [1 if cmodels.Permission.get_by_pk(x) in permission.deny.all() else 0 for x in range(1, 9)]]
        user_permissions_allow.append(user_perm_allow)
        user_permissions_deny.append(user_perm_deny)
      
    for permission in group_permissions_db:
        group_perm_allow = [permission.group.id, [1 if cmodels.Permission.get_by_pk(x) in permission.allow.all() else 0 for x in range(1, 9)]]       
        group_perm_deny = [permission.group.id, [1 if cmodels.Permission.get_by_pk(x) in permission.deny.all() else 0 for x in range(1, 9)]] 
        group_permissions_allow.append(group_perm_allow)
        group_permissions_deny.append(group_perm_deny)
    
    queryset = chelpers.objectToQuerySet(item) if item else None

    return chelpers.render_queryset_to_response(
                request,
                species = species,
                models = [model],
                queryset = queryset,
                template = "cyano/permission_edit.html" if edit else "cyano/permission.html",
                data = {
                    'users': users,
                    'groups': groups,
                    'user_permissions_allow': user_permissions_allow,
                    'group_permissions_allow': group_permissions_allow,
                    'user_permissions_deny': user_permissions_deny,
                    'group_permissions_deny': group_permissions_deny,
                    'permission_types': cmodels.Permission.permission_types
                })

@login_required
@resolve_to_objects
def jobs(request, species = None):
    from djcelery.models import TaskMeta
    from celery.task.control import inspect
    from ast import literal_eval # safer than eval
    
    pending = []
    finished = []
    running = []
    # Fetch pending tasks
    insp = inspect()
    message_debug = "No worker available"
    try:
        res = insp.reserved()
    except Exception as e:
        res = None
        message_debug = str(e)

    # Admins see all jobs
    is_admin = request.user.profile.is_admin()

    if not res:
        return chelpers.render_queryset_to_response_error(
            request,
            error = 503,
            msg = "Server error or no worker available. Please report this to an administrator!",
            msg_debug = message_debug
        )

    for v in res.values():
        for job in v:
            kwargs = literal_eval(job["kwargs"])
            if is_admin or kwargs["user"] == request.user.pk or kwargs["user"] == request.user.username:
                pending.append(kwargs)

    obj = TaskMeta.objects.all().order_by("pk")
    for o in obj:
        if is_admin or o.result["user"] == request.user.pk or o.result["user"] == request.user.username:
            if o.status == "SUCCESS" or o.status == "FAILURE":
                finished.append(o)
            elif o.status == "PROGRESS":
                running.append(o.result)

    finished = sorted(finished, key = lambda f: f.date_done, reverse = True)

    return chelpers.render_queryset_to_response(
        species = species,
        request = request, 
        models = [cmodels.UserProfile],
        template = 'cyano/jobs.html',
            data = {'pending': pending,
                    'finished': finished,
                    'running': running})

def sbgn(request):
    return chelpers.render_queryset_to_response(
        request,
        template="cyano/sbgn.html",
    )
