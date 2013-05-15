'''
Whole-cell knowledge base views

Affiliation: Covert Lab, Department of Bioengineering, Stanford University
Last updated: 2012-07-17
'''

from copy import deepcopy
from django.contrib.auth import login as auth_login, logout as auth_logout
from django.contrib.auth.decorators import login_required
from django.contrib.auth.forms import AuthenticationForm
from django.contrib.auth.models import User
from django.core.exceptions import ValidationError
from django.core.urlresolvers import reverse
from django.db.models import Count, Sum, Avg
from django.db.models.fields import BooleanField, NullBooleanField, AutoField, BigIntegerField, DecimalField, FloatField, IntegerField, PositiveIntegerField, PositiveSmallIntegerField, SmallIntegerField
from django.db.models.fields.related import OneToOneField, RelatedObject, ManyToManyField, ForeignKey
from django.db.models.query import EmptyQuerySet
from django.http import Http404, HttpResponse, HttpResponseRedirect
from django.shortcuts import get_object_or_404
from django.utils.text import capfirst
from django.views.decorators.debug import sensitive_post_parameters
from django.views.decorators.cache import never_cache
from django.views.decorators.csrf import csrf_protect
from haystack.query import SearchQuerySet
from itertools import chain
from cyano.forms import ExportDataForm, ImportDataForm
from cyano.helpers import getEntry, format_field_detail_view, objectToQuerySet, render_queryset_to_response, getObjectTypes, getModel, get_invalid_objects, get_edit_form_fields, get_edit_form_data
from cyano.helpers import validate_object_fields, validate_model_objects, validate_model_unique, save_object_data, batch_import_from_excel, readFasta
import cyano.models as models
from urlparse import urlparse
import numpy
import os
import settings
import tempfile

from cyano.helpers import render_queryset_to_response,\
    get_verbose_name_for_field_by_name, render_queryset_to_response_error

from cyano.decorators import resolve_to_objects
from cyano.decorators import permission_required

import cyano.models as models
from cyano.models import PermissionEnum as perm

def index(request):
    return render_queryset_to_response(
        request,
        template = "cyano/index.html",
        )

@resolve_to_objects
@permission_required(perm.READ_NORMAL)
def species(request, species):
    content = []
    if species is not None:        
        content.append([
            [0, 'Compartments', models.Compartment.objects.filter(species__id = species.id).count(), None, reverse('cyano.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Compartment'})],
        ])
        
        chrs = models.Chromosome.objects.filter(species__id = species.id)
        chrcontent = chrs.aggregate(length=Sum('length'));
        gc_content = 0 if len(chrs) == 0 else sum([chr.get_gc_content() * chr.length for chr in chrs]) / chrcontent['length']        
        content.append([
            [0, 'Chromosomes', chrs.count(), None, reverse('cyano.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Chromosome'})],
            [1, 'Length', chrcontent['length'], 'nt'],
            [1, 'GC-content', ('%0.1f' % (gc_content * 100)), '%'],
        ])
                
        tus = models.TranscriptionUnit.objects.filter(species__id = species.id).annotate(num_genes = Count('genes'))
        mons = tus.filter(num_genes__lte = 1)
        nPolys = tus.filter(num_genes__gt = 1).count()
        content.append([
            [0, 'Transcription units', tus.count(), None, reverse('cyano.views.list', kwargs={'species_wid': species.wid, 'model_type': 'TranscriptionUnit'})],            
            [1, 'Monocistrons', tus.count() - nPolys],
            [1, 'Polycistrons', nPolys],
        ])
        
        content.append([
            [0, 'Genes', models.Gene.objects.filter(species__id = species.id).count(), None, reverse('cyano.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Gene'})],
            [1, 'mRNA', models.Gene.objects.filter(species__id = species.id, type__wid='mRNA').count(), None, reverse('cyano.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Gene'}) + '?type=mRNA'],
            [1, 'rRNA', models.Gene.objects.filter(species__id = species.id, type__wid='rRNA').count(), None, reverse('cyano.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Gene'}) + '?type=rRNA'],
            [1, 'sRNA', models.Gene.objects.filter(species__id = species.id, type__wid='sRNA').count(), None, reverse('cyano.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Gene'}) + '?type=sRNA'],
            [1, 'tRNA', models.Gene.objects.filter(species__id = species.id, type__wid='tRNA').count(), None, reverse('cyano.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Gene'}) + '?type=tRNA'],
        ])
        
        content.append([
            [0, 'Chromosome features', 
                models.ChromosomeFeature.objects.filter(species__id = species.id).count(),
                None,
                reverse('cyano.views.list', kwargs={'species_wid': species.wid, 'model_type': 'ChromosomeFeature'}),
                ],
            [1, 'DnaA boxes', 
                models.ChromosomeFeature.objects.filter(species__id = species.id, type__parent__wid='ChromosomeFeature-DnaA_box').count(),
                ],                
            [1, 'Short tandem repeats', 
                models.ChromosomeFeature.objects.filter(species__id = species.id, type__parent__wid='ChromosomeFeature-Short_Tandem_Repeat').count(),
                ],
            [1, 'Other', 
                models.ChromosomeFeature.objects.filter(species__id = species.id).exclude(type__parent__wid='ChromosomeFeature-DnaA_box').exclude(type__parent__wid='ChromosomeFeature-Short_Tandem_Repeat').count()],
        ])
        
        content.append([
            [0, 'Metabolites', models.Metabolite.objects.filter(species__id = species.id).count(), None, reverse('cyano.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Metabolite'})],
            [1, 'Amino acids', 
                models.Metabolite.objects.filter(species__id = species.id, type__wid='amino_acid').count() + 
                models.Metabolite.objects.filter(species__id = species.id, type__wid='modified_amino_acid').count() +
                models.Metabolite.objects.filter(species__id = species.id, type__wid='non-standard_amino_acid').count() +
                models.Metabolite.objects.filter(species__id = species.id, type__wid='vitamin_non-standard_amino_acid').count()                
                ],
            [1, 'Antibiotic', 
                models.Metabolite.objects.filter(species__id = species.id, type__wid='antibiotic').count() + 
                models.Metabolite.objects.filter(species__id = species.id, type__parent__wid='antibiotic').count()
                ],
            [1, 'Gases', 
                models.Metabolite.objects.filter(species__id = species.id, type__wid='gas').count() + 
                models.Metabolite.objects.filter(species__id = species.id, type__parent__wid='gas').count()
                ],
            [1, 'Ions', 
                models.Metabolite.objects.filter(species__id = species.id, type__wid='ion').count() + 
                models.Metabolite.objects.filter(species__id = species.id, type__parent__wid='ion').count()
                ],
            [1, 'Lipids', 
                models.Metabolite.objects.filter(species__id = species.id, type__wid='lipid').count() + 
                models.Metabolite.objects.filter(species__id = species.id, type__parent__wid='lipid').count() +
                models.Metabolite.objects.filter(species__id = species.id, type__parent__parent__wid='lipid').count()
                ],
            [1, 'Vitamins', 
                models.Metabolite.objects.filter(species__id = species.id, type__wid='vitamin').count() + 
                models.Metabolite.objects.filter(species__id = species.id, type__parent__wid='vitamin').count()
                ],
        ])
                
        mons = models.ProteinMonomer.objects.filter(species__id = species.id)
        cpxs = models.ProteinComplex.objects.filter(species__id = species.id)
        monDNABind = mons.filter(dna_footprint__length__gt=0).count()
        monIntMem = mons.filter(localization__wid = 'm').exclude(signal_sequence__type = 'lipoprotein').count() + mons.filter(localization__wid = 'tm').exclude(signal_sequence__type = 'lipoprotein').count()
        monLipo = mons.filter(signal_sequence__type = 'lipoprotein').count()
        monSecreted = mons.filter(signal_sequence__type = 'secretory').count()
        monTermOrg = mons.filter(localization__wid = 'tc').count() + mons.filter(localization__wid = 'tm').count()
        cpxDNABind = cpxs.filter(dna_footprint__length__gt=0).count()
        content.append([
            [0, 'Proteins', mons.count() + cpxs.count()],
                [1, 'Monomers', mons.count(), None, reverse('cyano.views.list', kwargs={'species_wid': species.wid, 'model_type': 'ProteinMonomer'})],            
                    [2, 'DNA-binding', monDNABind],
                    [2, 'Integral membrane', monIntMem],
                    [2, 'Lipoprotein', monLipo, None, reverse('cyano.views.list', kwargs={'species_wid': species.wid, 'model_type': 'ProteinMonomer'}) + '?signal_sequence__type=lipoprotein'],            
                    [2, 'Secreted', monSecreted, None, reverse('cyano.views.list', kwargs={'species_wid': species.wid, 'model_type': 'ProteinMonomer'}) + '?signal_sequence__type=secretory'],
                    [2, 'Terminal organelle', monTermOrg],                    
                [1, 'Complexes', cpxs.count(), None, reverse('cyano.views.list', kwargs={'species_wid': species.wid, 'model_type': 'ProteinComplex'})],
                    [2, 'DNA-binding', cpxDNABind],
        ])

        rxns = models.Reaction.objects.filter(species__id = species.id)
        content.append([
            [0, 'Reactions', rxns.count(), None, reverse('cyano.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'})],
            [1, 'DNA damage', rxns.filter(processes__wid='Process_DNADamage').count(), None, reverse('cyano.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_DNADamage'],
            [1, 'DNA repair', rxns.filter(processes__wid='Process_DNARepair').count(), None, reverse('cyano.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_DNARepair'],
            [1, 'Metabolic', rxns.filter(processes__wid='Process_Metabolism').count(), None, reverse('cyano.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_Metabolism'],            
            [1, 'Protein decay', rxns.filter(processes__wid='Process_ProteinDecay').count(), None, reverse('cyano.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_ProteinDecay'],
            [1, 'Protein modification', rxns.filter(processes__wid='Process_ProteinModification').count(), None, reverse('cyano.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_ProteinModification'],            
            [1, 'Replication Initiation', rxns.filter(processes__wid='Process_ReplicationInitiation').count(), None, reverse('cyano.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_ReplicationInitiation'],
            [1, 'RNA decay', rxns.filter(processes__wid='Process_RNADecay').count(), None, reverse('cyano.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_RNADecay'],
            [1, 'RNA modification', rxns.filter(processes__wid='Process_RNAModification').count(), None, reverse('cyano.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_RNAModification'],            
            [1, 'RNA processing', rxns.filter(processes__wid='Process_RNAProcessing').count(), None, reverse('cyano.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_RNAProcessing'],            
            [1, 'Transcription', rxns.filter(processes__wid='Process_Transcription').count(), None, reverse('cyano.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_Transcription'],            
            [1, 'Translation', rxns.filter(processes__wid='Process_Translation').count(), None, reverse('cyano.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_Translation'],            
            [1, 'tRNA aminoacylation', rxns.filter(processes__wid='Process_tRNAAminoacylation').count(), None, reverse('cyano.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_tRNAAminoacylation'],
            [1, 'Other', rxns.exclude(processes__wid='Process_DNADamage')
                .exclude(processes__wid='Process_DNARepair')
                .exclude(processes__wid='Process_Metabolism')
                .exclude(processes__wid='Process_ProteinDecay')
                .exclude(processes__wid='Process_ProteinModification')
                .exclude(processes__wid='Process_ReplicationInitiation')                
                .exclude(processes__wid='Process_RNADecay')
                .exclude(processes__wid='Process_RNAModification')
                .exclude(processes__wid='Process_RNAProcessing')
                .exclude(processes__wid='Process_Transcription')
                .exclude(processes__wid='Process_Translation')
                .exclude(processes__wid='Process_tRNAAminoacylation')
                .count()],
        ])        
        
        tr = models.TranscriptionalRegulation.objects.filter(species__id = species.id)
        nTus = len(set([x[0] for x in tr.values_list('transcription_unit')]))
        nTfs = len(set([x[0] for x in tr.values_list('transcription_factor')]))
        content.append([
            [0, 'Transcriptional regulation'],
            [1, 'Interactions', tr.count(), None, reverse('cyano.views.list', kwargs={'species_wid': species.wid, 'model_type': 'TranscriptionalRegulation'})],
            [1, 'Transcriptional regulators', nTfs],
            [1, 'Regulated promoters', nTus],
        ])        

        content.append([
            [0, 'Pathways', models.Pathway.objects.filter(species__id = species.id).count(), None, reverse('cyano.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Pathway'})],
        ])
        content.append([
            [0, 'Stimuli', models.Stimulus.objects.filter(species__id = species.id).count(), None, reverse('cyano.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Stimulus'})],
        ])
        
        
        nCellComp = models.Metabolite.objects.filter(species__id = species.id, biomass_composition__concentration__isnull=False).count()
        nMediaComp = models.Metabolite.objects.filter(species__id = species.id, media_composition__concentration__isnull=False).count()        
        nKineticsKeq = models.Reaction.objects.filter(species__id = species.id, keq__isnull=False).count()
        nKineticsKm = \
            models.Reaction.objects.filter(species__id = species.id, kinetics_forward__km__isnull=False).count() + \
            models.Reaction.objects.filter(species__id = species.id, kinetics_backward__km__isnull=False).count()
        nKineticsVmax = \
            models.Reaction.objects.filter(species__id = species.id, kinetics_forward__vmax__isnull=False).count() + \
            models.Reaction.objects.filter(species__id = species.id, kinetics_backward__vmax__isnull=False).count()
        nRnaExp = models.Gene.objects.filter(species__id = species.id, expression__isnull=False).count()
        nRnaHl = models.Gene.objects.filter(species__id = species.id, half_life__isnull=False).count()
        nStimuli = models.Stimulus.objects.filter(species__id = species.id, value__isnull=False).count()
        nTrAffinity = models.TranscriptionalRegulation.objects.filter(species__id = species.id, affinity__isnull=False).count()        
        nTrActivity = models.TranscriptionalRegulation.objects.filter(species__id = species.id, activity__isnull=False).count()        
        nOther = models.Parameter.objects.filter(species__id = species.id).count()
        nTotParameters = nCellComp + nMediaComp + nKineticsVmax + nRnaExp + nRnaHl + nStimuli + nTrAffinity + nTrActivity + nOther

        content.append([
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
            [1, 'Other', nOther, None, reverse('cyano.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Parameter'})],            
        ])
        
        content.append([
            [0, 'Processes', models.Process.objects.filter(species__id = species.id).count(), None, reverse('cyano.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Process'})],
        ])
        
        content.append([
            [0, 'States', models.State.objects.filter(species__id = species.id).count(), None, reverse('cyano.views.list', kwargs={'species_wid': species.wid, 'model_type': 'State'})],
        ])
        
    nContent = [len(x) for x in content]
    totContent = sum(nContent)
    cum = 0
    idx = 0
    breakIdxs = [0, 0]
    for x in nContent:
        cum += x
        idx += 1
        if cum > totContent * 1/ 3 and breakIdxs[0] == 0:
            breakIdxs[0] = idx
        if cum > totContent * 2 / 3 and breakIdxs[1] == 0:
            breakIdxs[1] = idx        
            
    contentCol1 = []
    contentCol2 = []
    contentCol3 = []
    i = 0
    for x in content[:breakIdxs[0]]:
        i += 1
        for y in x:
            contentCol1.append([i] + y)    
    i = 0
    for x in content[breakIdxs[0]:breakIdxs[1]]:
        i += 1
        for y in x:
            contentCol2.append([i] + y)
    i = 0
    for x in content[breakIdxs[1]:]:
        i += 1
        for y in x:
            contentCol3.append([i] + y)
        
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
    if species is not None:
        #refs = models.PublicationReference.objects.filter(species__id = species.id)
        #sources['total'] = refs.count()
        #sources['types'] = [
        #        {'type': 'Articles', 'count': refs.filter(species__id = species.id, type__wid='article').count()},
        #        {'type': 'Books', 'count': refs.filter(species__id = species.id, type__wid='book').count()},
        #        {'type': 'Thesis', 'count': refs.filter(species__id = species.id, type__wid='thesis').count()},
        #        {'type': 'Other', 'count': refs.filter(species__id = species.id, type__wid='misc').count()},
        #    ]
        #sources['dates'] = refs.filter(year__isnull=False).order_by('year').values('year').annotate(count=Count('year'))
        
            
        nEstimated = models.Parameter.objects.filter(species__id = species.id, value__evidence__is_experimentally_constrained=False).count()
        nExpConstrained = nTotParameters - nEstimated
        sources['evidence_parameters'] = [
            {'type': 'Experimentally constrained', 'count': nExpConstrained},
            {'type': 'Computationally estimated', 'count': nEstimated},
            ]
            
        
        sources['evidence_species'] = models.Evidence.objects.filter(species_component__species__id = species.id).values('species').annotate(count = Count('id'))
        sources['evidence_media'] = models.Evidence.objects.filter(species_component__species__id = species.id).values('media').annotate(count = Count('id'))
        sources['evidence_pH'] = models.Evidence.objects.filter(species_component__species__id = species.id).values('pH').annotate(count = Count('id'))
        sources['evidence_temperature'] = models.Evidence.objects.filter(species_component__species__id = species.id).values('temperature').annotate(count = Count('id'))
            
    return render_queryset_to_response(
        species = species,
        data = {
            'content': [contentCol1, contentCol2, contentCol3],
            'contentRows': range(max(len(contentCol1), len(contentCol2), len(contentCol3))),
            'sources': sources,            
            },
        request = request, 
        template = 'cyano/species.html')
    
def about(request, species_wid=None):
    return render_queryset_to_response(
        species = species,
        request = request, 
        template = 'cyano/about.html', 
        data = {
            'ROOT_URL': settings.ROOT_URL,
        }
    )    
        
def tutorial(request, species_wid=None):
    return render_queryset_to_response(
        species = species,
        request = request, 
        template = 'cyano/tutorial.html')    

@login_required
@resolve_to_objects
def users(request, species = None):
    queryset = models.UserProfile.objects.all().filter(user__is_active = True)
    return render_queryset_to_response(
        species = species,
        request = request, 
        models = [models.UserProfile],
        queryset = queryset,
        template = 'cyano/users.html')    
        
@login_required   
@resolve_to_objects
def user(request, username, species = None):
    queryset = objectToQuerySet(get_object_or_404(models.UserProfile, user__username = username), model = models.UserProfile)
    return render_queryset_to_response(
        species = species,
        request = request,
        models = [models.UserProfile],
        queryset = queryset,
        template = 'cyano/user.html')    

def search(request, species_wid = None):
    query = request.GET.get('q', '')
    engine = request.GET.get('engine', 'haystack')
    
    if engine == 'haystack':
        return search_haystack(request, species_wid, query)
    else:
        return search_google(request, species_wid, query)
                
def search_haystack(request, species_wid, query):
    #search
    if species_wid is None:
        species_wid = models.Species.objects.all()[0].wid
    results = SearchQuerySet().filter(species=species).filter(content=query)
    
    #calculate facets        
    facets = results.facet('model_type')
    tmp = facets.facet_counts()['fields']['model_type']
    modelNameFacet = []
    objectTypes = getObjectTypes()
    models = []
    for tmp2 in tmp:
        modelName = objectTypes[objectTypes.index(tmp2[0])]
        modelNameFacet.append({
            'name':modelName, 
            'verbose_name': getModel(modelName)._meta.verbose_name,
            'count':tmp2[1],
            })
        models.append(getModel(modelName))
    modelNameFacet.sort(lambda x, y:cmp(x['verbose_name'], y['verbose_name']))
    
    #narrow search by facets
    model_type = request.GET.get('model_type', '')
    if model_type:
        results = results.models(getModel(model_type))
        
    #order results
    results = results.order_by('wid')
    
    #convert results to query set
    queryset = EmptyQuerySet()
    for object in results:
        tmp = object.model.objects.none()
        tmp._result_cache.append(object.object)
        queryset = chain(queryset, tmp)
    
    #form response
    return render_queryset_to_response(
        species = species,
        request = request, 
        models = models, 
        queryset = queryset, 
        template = 'cyano/search.html', 
        data = {
            'query': query,
            'engine': 'haystack',
            'model_type': model_type,
            'modelNameFacet': modelNameFacet,
            })

def search_google(request, species_wid, query):
    return render_queryset_to_response(
        species = species,
        request = request, 
        template = 'cyano/googleSearch.html', 
        data = {
            'query': query,
            'engine': 'google',
            })

@resolve_to_objects
#@permission_required(perm.READ_NORMAL)
def list(request, species, model):
    #try:
    #    getObjectTypes().index(model_type)
    #except ValueError:
    #    raise Http404

    objects = model.objects.all().filter(species__id=species.id)
    
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
                
        if isinstance(field, (ForeignKey, ManyToManyField)) and not issubclass(field.rel.to, models.Entry):
            continue
        
        if isinstance(field, (ForeignKey, ManyToManyField)):
            tmp = model.objects.filter(species__id=species.id).order_by(field_full_name + '__name').values(field_full_name).annotate(count=Count(field_full_name))
        else:
            tmp = model.objects.filter(species__id=species.id).order_by(field_full_name).values(field_full_name).annotate(count=Count(field_full_name))
        facets = []
        for facet in tmp:
            value = facet[field_full_name]            
            if value is None or unicode(value) == '':
                continue
            
            if isinstance(field, (ForeignKey, ManyToManyField)):
                tmp2 = tmp_model.objects.values('wid', 'name').get(id=value)
                id = tmp2['wid']
                name = capfirst(tmp2['name'])
            elif (field.choices is not None) and (len(field.choices) > 0) and (not isinstance(field, (BooleanField, NullBooleanField))):    
                id = value
                choices = [x[0] for x in field.choices]
                if id in choices:
                    name = field.choices[choices.index(id)][1]
                else:
                    name = capfirst(value)
            else:
                id = value
                name = capfirst(value)
            if value is not None and unicode(value) != '':
                facets.append({
                    'id': unicode(id), 
                    'name': unicode(name),
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

    return render_queryset_to_response(
        species = species,
        request = request, 
        models = [model], 
        queryset = objects, 
        template = 'cyano/list.html', 
        data = {
            'facet_fields': facet_fields,
            })

@resolve_to_objects
#@permission_required(perm.READ_NORMAL)
def detail(request, species, model, item):
    fieldsets = deepcopy(model._meta.fieldsets)
    
    #filter out type, metadata
    fieldset_names = [x[0] for x in fieldsets]
    if 'Type' in fieldset_names:
        idx = fieldset_names.index('Type')
        del fieldsets[idx]
        
    #filter out empty fields
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
                
            data = format_field_detail_view(species, item, field_name, request.user.is_anonymous())
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
    qs = objectToQuerySet(item, model = model)

    #render response
    return render_queryset_to_response(
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
#@permission_required(perm.READ_HISTORY)
def history(request, species, model = None, item = None, detail_id = None):
    if item:    
        objects = objectToQuerySet(item, model = model)[0].revisions.prefetch_related("detail", "detail__user", "current", "current__model_type").order_by("-detail__id").distinct("detail__id")
    elif model:
        objects = model.objects.filter(species__id=species.id)#.prefetch_related("revisions", "revisions__detail", "revisions__detail__user", "model_type").order_by("-revisions__detail__id").distinct("revisions__detail__id")
        #tmeta = models.TableMeta.objects.get(name = model._meta.object_name)
        #objects = models.Revision.objects.filter(current__)
    else:
        objects = models.SpeciesComponent.objects.filter(species__id=species.id).prefetch_related("detail", "detail__user", "current", "current__model_type").order_by("-detail__id").distinct("detail__id")

    #objects = objects.prefetch_related("detail", "detail__user", "current", "current__model_type").order_by("-detail__id").distinct("detail__id")

    revisions = []
    
    entry = []
    date = None
    
    for revision in objects:
        last_date = date
        wid = revision.current.wid
        item_model = revision.current.model_type.name
        detail_id = revision.detail.pk
        date = revision.detail.date.date()
        time = revision.detail.date.strftime("%H:%M")
        reason = revision.detail.reason
        author = revision.detail.user
        url = reverse("cyano.views.history", kwargs = {"species_wid": species.wid, "model_type": item_model, "wid": wid, "detail_id": detail_id})
        
        if last_date != date:
            revisions.append(entry)
            entry = [date, []]
        
        entry[1].append({'id': detail_id, 'time': time, 'wid': wid, 'reason': reason, 'author': author, 'url': url})
    revisions.append(entry)
    
    if item:
        qs = objectToQuerySet(item, model = model)
    else:
        qs = objects

    return render_queryset_to_response(
        species = species,
        request = request,
        models = [model],
        queryset = qs,
        template = 'cyano/history.html',
        data = {
            'revisions': revisions[1:]
            })
            
@login_required        
def add(request, model_type, species_wid=None):
    return edit(request, model_type=model_type, species=species, action='add')
        
@login_required
def edit(request, wid=None, model_type=None, species_wid=None, action='edit'):
    #retrieve object
    if action == 'edit':
        obj = getEntry(species = species, wid = wid)
        if obj is None:
            raise Http404
        model = obj.__class__ 
    else:        
        model = getModel(model_type)
        obj = model()
    
    #save object
    error_messages = {}
    if request.method == 'POST':
        submitted_data = get_edit_form_data(model, request.POST)
        
        data = submitted_data
        data['id'] = obj.id
        data['species'] = species_wid
        data['model_type'] = model.__name__
        
        try:
            #validate is WID unique
            if issubclass(model, models.SpeciesComponent):
                qs = models.SpeciesComponent.objects.values('wid', 'model_type').filter(species__wid=species_wid)
            else:
                qs = model.objects.values('wid', 'model_type').all()
                
            if action == 'edit':
                qs = qs.exclude(id=obj.id)
                
            wids = {}
            for x in qs:
                wids[x['wid']] = x['model_type']
            
            if data['wid'] in wids.keys():
                raise ValidationError({'wid': 'Value must be unique'})
                
            wids[data['wid']] = model.__name__
        
            #validate
            data = validate_object_fields(model, data, wids, species_wid, data['wid'])
            validate_model_objects(model, data)
            validate_model_unique(model, [data])
            
            #save
            obj = save_object_data(species_wid, obj, data, {}, request.user, save=False, save_m2m=False)
            obj = save_object_data(species_wid, obj, data, {data['wid']: obj}, request.user, save=True, save_m2m=False)
            obj = save_object_data(species_wid, obj, data, {data['wid']: obj}, request.user, save=True, save_m2m=True)
            
            #redirect to details page
            return HttpResponseRedirect(obj.get_absolute_url())
        except ValidationError as error:
            error_messages = error.message_dict
    
    #form query set
    if action == 'edit':
        obj = getEntry(species = species, wid = wid)
        if obj is None:
            raise Http404
        qs = objectToQuerySet(obj, model = model)
    else:
        obj = None
        qs = model.objects.none()
        
    #display form
    fields, initial_values = get_edit_form_fields(species_wid, model, obj=obj)
    
    if request.method == 'POST':
        initial_values = submitted_data
    return render_queryset_to_response(
        species = species,
        request = request, 
        models = [model],
        queryset = qs,
        template = 'cyano/edit.html', 
        data = {
            'model_verbose_name': model._meta.verbose_name,
            'action': action,
            'fields': fields,
            'references_choices': models.PublicationReference.objects.filter(species__wid = species_wid).values_list('wid'),
            'initial_values': initial_values,
            'error_messages': error_messages,
            }
        )
            
@login_required        
def delete(request, species_wid, wid):
    #retrieve object
    obj = getEntry(species = species, wid = wid)
    if obj is None:
        raise Http404    
    model = obj.__class__ 
    qs = objectToQuerySet(obj, model = model)
    
    #delete
    if request.method == 'POST':
        obj.delete()
        return HttpResponseRedirect(reverse('cyano.views.list', kwargs={'species_wid':species_wid, 'model_type': model.__name__}))
        
    #confirmation message
    return render_queryset_to_response(
        species = species,
        request = request, 
        models = [model],
        queryset = qs,
        template = 'cyano/delete.html', 
        data = {
            'model_verbose_name': model._meta.verbose_name
            }
        )

def exportData(request, species_wid=None):
    getDict = request.GET.copy()
    if getDict.get('format', ''):
        getDict.__setitem__('species', getDict.get('species', species_wid))    
    form = ExportDataForm(getDict or None)
    if not form.is_valid():        
        return render_queryset_to_response(
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
            model_types = getObjectTypes()
        else:
            model_types = form.cleaned_data['model_type']
        
        for model_type in model_types:
            model = getModel(model_type)
            if issubclass(model, models.SpeciesComponent):
                queryset = chain(queryset, model.objects.filter(species__id=species.id).select_related(depth=2).all())
            else:
                queryset = chain(queryset, model.objects.select_related(depth=2).filter(id=species.id))
            models.append(getModel(model_type))
        
        return render_queryset_to_response(
            species = species,
            request = request, 
            queryset = queryset, 
            template = 'cyano/exportDataResult.html', 
            models = models)

@login_required
@resolve_to_objects
def importData(request, species=None):
    import cyano.importer.fasta as FastaImporter
    from cyano.importer.genbank import Genbank as GenbankImporter
    if request.method == 'POST':
        form = ImportDataForm(request.POST or None, request.FILES)
        
        if form.is_valid():       
            selected_species_wid = form.cleaned_data['species'] or form.cleaned_data['new_wid']
            if selected_species_wid == '' or selected_species_wid is None:
                form.errors["species"] = ['Please select a species or create a new one']
            elif form.cleaned_data['new_wid'] and models.Species.objects.filter(wid = selected_species_wid).exists():
                form.errors["new_wid"] = ['The identifier specified is already in use']
            else:
                #save to temporary file
                originalFileName, originalFileExtension = os.path.splitext(request.FILES['file'].name)
                fid = tempfile.NamedTemporaryFile(suffix = originalFileExtension, delete = False)
                filename = fid.name
                for chunk in request.FILES['file'].chunks():
                    fid.write(chunk)
                fid.close()
                
                #read file
                data_type = form.cleaned_data["data_type"]
                request.POST["new_species"] = form.cleaned_data["new_species"]
                request.POST["data_type"] = data_type
                request.POST["reason"] = form.cleaned_data["reason"]
                
                if data_type == "fasta":
                    f = FastaImporter.Fasta()
                    f.load(filename)
                    return f.preview(request, selected_species_wid)
                elif data_type == "fastagene":
                    g = GenbankImporter()
                    g.load(filename)
                    return g.preview(request, selected_species_wid)
                
                os.remove(filename)
                raise Http404
                
                #if originalFileExtension == '.xlsx':
                #    try:
                #        batch_import_from_excel(selected_species_wid, filename, request.user)
                #        success = True
                #        message = 'Data successfully saved!'
                #    except ValidationError as error:
                #        success = False
                #        message = 'Unable to import data: ' + ' '.join(error.messages)

    else:
        form = ImportDataForm(None)

    return render_queryset_to_response(
        species = species,
        request = request, 
        template = 'cyano/importDataForm.html', 
        data = {
            'form': form},
        )

@login_required
@resolve_to_objects
def importSubmitData(request, species=None):
    from cyano.importer.genbank import Genbank as GenbankImporter
    if request.method == 'POST':
        data_type = request.POST["data_type"]
        filename = request.POST['filename']
        #if data_type == "fasta":
            
        if data_type == "fastagene":
            g = GenbankImporter()
            g.load(filename)
            os.remove(filename)
            return g.submit(request, request.POST['species_wid'])

        os.remove(filename)
        raise Http404
    pass
    
def validate(request, species_wid):
    errors = get_invalid_objects(models.Species.objects.values('id').get(wid=species_wid)['id'])
    
    return render_queryset_to_response(
        species = species,
        request = request, 
        template = 'cyano/validate.html', 
        data = {            
            'errors': errors
            },
        )
        
@sensitive_post_parameters()
@csrf_protect
@never_cache
@resolve_to_objects
def login(request, species=None):
    next = request.REQUEST.get('next', '')
    
    if request.method == "POST":
        form = AuthenticationForm(data=request.POST)
        if form.is_valid():
            auth_login(request, form.get_user())
            
            if request.session.test_cookie_worked():
                request.session.delete_test_cookie()

            return HttpResponseRedirect(next)
    else:
        form = AuthenticationForm(request)

    request.session.set_test_cookie()

    return render_queryset_to_response(
        species = species,
        request = request, 
        template = 'cyano/login.html', 
        data = {
            'form': form,
            'next': next,
        })

@resolve_to_objects        
def logout(request, species=None):
    auth_logout(request)    
    return render_queryset_to_response(
        species = species,
        request = request, 
        template = 'cyano/logout.html', 
        )
    
def sitemap(request):
    return render_queryset_to_response(
        request = request, 
        template = 'cyano/sitemap.xml', 
        data = {
            'ROOT_URL': settings.ROOT_URL,
            'qs_species': models.Species.objects.all(),
        }
    )
    
def sitemap_toplevel(request):
    return render_queryset_to_response(
        request = request, 
        template = 'cyano/sitemap_toplevel.xml', 
        data = {
            'ROOT_URL': settings.ROOT_URL,
        }
    )

@resolve_to_objects
def sitemap_species(request, species):
    return render_queryset_to_response(
        request = request, 
        template = 'cyano/sitemap_species.xml', 
        data = {
            'ROOT_URL': settings.ROOT_URL,
            'entries': models.SpeciesComponent.objects.filter(species__id = species.id),
        }
    )

@resolve_to_objects
#@permission_required(perm.READ_PERMISSION)
def permission(request, species, model = None, item = None):
    users = models.UserProfile.objects.all().filter(user__is_active = True)
    groups = models.GroupProfile.objects.all()
    
    # Permissions for item or species depending on page
    entry = species if item == None else item

    if request.method == 'POST':
        # Form submit -> Save changes
        post_types = models.Permission.permission_types + ["uid_allow", "gid_allow", "uid_deny", "gid_deny"]
        error = False
        
        # Check that all needed fields are in the request
        if not all(p in request.POST.keys() for p in post_types):
            error = True

        if not error:
            new_permissions = {}

            # Remove all fields except the needed ones, verify
            # that all new permissions are numbers
            for k, v in request.POST.lists():
                if k in models.Permission.permission_types + post_types:
                    try:
                        new_permissions[k] = map(lambda x: int(x), v)
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
            user_verify = models.UserProfile.objects.filter(pk__in = new_permissions["uid_allow"]).count()
            user_verify += models.UserProfile.objects.filter(pk__in = new_permissions["uid_deny"]).count()
            group_verify = models.GroupProfile.objects.filter(pk__in = new_permissions["gid_allow"]).count()
            group_verify += models.GroupProfile.objects.filter(pk__in = new_permissions["gid_deny"]).count()

            # If number of submitted data larger then database data at least one item was invalid
            if user_count_allow + group_count_allow + user_count_deny + group_count_deny != user_verify + group_verify:
                error = True

        if error:
            return render_queryset_to_response_error(
                request,
                species = species,
                error = 400,
                msg = "Invalid POST request on permission page.",
                msg_debug = "POST was: " + str(request.POST.lists()))        
        
        # Data is valid, begin with database operations
        
        new_permissions_user_allow = {}
        new_permissions_group_allow = {}
        new_permissions_user_deny = {}
        new_permissions_group_deny = {}

        # Rotate array to user -> permissions
        for i, user in enumerate(new_permissions["uid_allow"]):
            new_permissions_user_allow[models.UserProfile.objects.get(pk = user)] = [new_permissions[x][i] for x in models.Permission.permission_types]
        
        # Rotate array to group -> permissions
        for i, group in enumerate(new_permissions["gid_allow"], user_count_allow):  
            new_permissions_group_allow[models.GroupProfile.objects.get(pk = group)] = [new_permissions[x][i] for x in models.Permission.permission_types]
        
        # Rotate array to user -> permissions
        for i, user in enumerate(new_permissions["uid_deny"], user_count_allow + group_count_allow):
            new_permissions_user_deny[models.UserProfile.objects.get(pk = user)] = [new_permissions[x][i] for x in models.Permission.permission_types]
        
        # Rotate array to group -> permissions
        for i, group in enumerate(new_permissions["gid_deny"], user_count_allow + group_count_allow + user_count_deny):  
            new_permissions_group_deny[models.GroupProfile.objects.get(pk = group)] = [new_permissions[x][i] for x in models.Permission.permission_types]
        
        # Check if user or group is completly missing but had permissions before
        # Delete if this is the case
        all_user_perms = models.UserPermission.objects.filter(entry = entry).prefetch_related("user")
        for permission in all_user_perms:
            if not permission.user in new_permissions_user_allow and not permission.user in new_permissions_user_deny:
                permission.delete()
            else:
                if not permission.user in new_permissions_user_allow:
                    permission.allow.clear()
                if not permission.user in new_permissions_user_deny:
                    permission.deny.clear()
        
        all_group_perms = models.GroupPermission.objects.filter(entry = entry).prefetch_related("group")
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
            new_perm, _ = models.UserPermission.objects.get_or_create(entry = entry, user = u)
            allows = new_perm.allow.all()
            for i, x in enumerate(range(1, 9)):
                if p[i] == 0:
                    # delete the entry if it exists
                    if models.Permission.get_by_pk(x) in allows:
                        new_perm.allow.remove(models.Permission.get_by_pk(x))
                else:
                    # add the entry if missing
                    if not models.Permission.get_by_pk(x) in allows:
                        new_perm.allow.add(models.Permission.get_by_pk(x))

        for u, p in new_permissions_user_deny.iteritems():
            new_perm, _ = models.UserPermission.objects.get_or_create(entry = entry, user = u)
            denies = new_perm.deny.all()
            for i, x in enumerate(range(1, 9)):
                if p[i] == 0:
                    # delete the entry if it exists
                    if models.Permission.get_by_pk(x) in denies:
                        new_perm.deny.remove(models.Permission.get_by_pk(x))
                else:
                    # add the entry if missing
                    if not models.Permission.get_by_pk(x) in denies:
                        new_perm.deny.add(models.Permission.get_by_pk(x))

            new_perm.save()
                    
        for g, p in new_permissions_group_allow.iteritems():
            new_perm, _ = models.GroupPermission.objects.get_or_create(entry = entry, group = g)
            allows = new_perm.allow.all()
            for i, x in enumerate(range(1, 9)):
                if p[i] == 0:
                    # delete the entry if it exists
                    if models.Permission.get_by_pk(x) in allows:
                        new_perm.allow.remove(models.Permission.get_by_pk(x))
                else:
                    # add the entry if missing
                    if not models.Permission.get_by_pk(x) in allows:
                        new_perm.allow.add(models.Permission.get_by_pk(x))
            
        for g, p in new_permissions_group_deny.iteritems():     
            new_perm, _ = models.GroupPermission.objects.get_or_create(entry = entry, group = g)    
            denies = new_perm.deny.all()
            for i, x in enumerate(range(1, 9)):
                if p[i] == 0:
                    # delete the entry if it exists
                    if models.Permission.get_by_pk(x) in denies:
                        new_perm.deny.remove(models.Permission.get_by_pk(x))
                else:
                    # add the entry if missing
                    if not models.Permission.get_by_pk(x) in denies:
                        new_perm.deny.add(models.Permission.get_by_pk(x))
            
            new_perm.save()

    # Rendering of the permission page
    user_permissions_db = models.UserPermission.objects.filter(entry = entry).prefetch_related("allow", "deny", "user", "user__user")
    group_permissions_db = models.GroupPermission.objects.filter(entry = entry).prefetch_related("allow", "deny", "group", "group__group")

    user_permissions_allow = []
    group_permissions_allow = []
    user_permissions_deny = []
    group_permissions_deny = []
    
    # Calculate permission array for user and groups
    # Contains user mapping to 8 numbers, 1 if permission available, 0 if not
    for permission in user_permissions_db:
        user_perm_allow = [permission.user.id, [1 if models.Permission.get_by_pk(x) in permission.allow.all() else 0 for x in range(1, 9)]]
        user_perm_deny = [permission.user.id, [1 if models.Permission.get_by_pk(x) in permission.deny.all() else 0 for x in range(1, 9)]]
        user_permissions_allow.append(user_perm_allow)
        user_permissions_deny.append(user_perm_deny)
      
    for permission in group_permissions_db:
        group_perm_allow = [permission.group.id, [1 if models.Permission.get_by_pk(x) in permission.allow.all() else 0 for x in range(1, 9)]]       
        group_perm_deny = [permission.group.id, [1 if models.Permission.get_by_pk(x) in permission.deny.all() else 0 for x in range(1, 9)]] 
        group_permissions_allow.append(group_perm_allow)
        group_permissions_deny.append(group_perm_deny)

    return render_queryset_to_response(
                request,
                species = species,
                template = "cyano/permission.html",
                data = {
                    'users': users,
                    'groups': groups,
                    'user_permissions_allow': user_permissions_allow,
                    'group_permissions_allow': group_permissions_allow,
                    'user_permissions_deny': user_permissions_deny,
                    'group_permissions_deny': group_permissions_deny,
                    'permission_types': models.Permission.permission_types
                })