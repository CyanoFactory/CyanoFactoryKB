'''
Whole-cell knowledge base views

Affiliation: Covert Lab, Department of Bioengineering, Stanford University
Last updated: 2012-07-17
'''

import os
import settings
import tempfile
from copy import deepcopy
from itertools import chain

from django.contrib.auth import login as auth_login, logout as auth_logout
from django.contrib.auth.decorators import login_required
from django.contrib.auth.forms import AuthenticationForm
from django.core.exceptions import ValidationError
from django.core.urlresolvers import reverse
from django.db import transaction
from django.db.models import Count, Sum
from django.db.models.fields import BooleanField, NullBooleanField, AutoField, BigIntegerField, DecimalField, FloatField, IntegerField, PositiveIntegerField, PositiveSmallIntegerField, SmallIntegerField
from django.db.models.fields.related import RelatedObject, ManyToManyField, ForeignKey
from django.db.models.query import EmptyQuerySet
from django.http import Http404, HttpResponseRedirect
from django.shortcuts import get_object_or_404
from django.utils.text import capfirst
from django.views.decorators.debug import sensitive_post_parameters
from django.views.decorators.cache import never_cache
from django.views.decorators.csrf import csrf_protect

from haystack.query import SearchQuerySet

from cyano.forms import ExportDataForm, ImportDataForm
import cyano.helpers as chelpers
import cyano.models as cmodels
from cyano.models import PermissionEnum as perm
from cyano.decorators import resolve_to_objects, permission_required

def index(request):
    return chelpers.render_queryset_to_response(
        request,
        template = "cyano/index.html",
        )

@resolve_to_objects
@permission_required(perm.READ_NORMAL)
def species(request, species):
    content = []
    if species is not None:
        content.append([
            [0, 'Compartments', cmodels.Compartment.objects.filter(species__id = species.id).count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Compartment'})],
        ])
        
        chrs = cmodels.Chromosome.objects.filter(species__id = species.id)
        chrcontent = chrs.aggregate(length=Sum('length'));
        gc_content = 0 if len(chrs) == 0 else sum([chro.get_gc_content() * chro.length for chro in chrs]) / chrcontent['length']        
        content.append([
            [0, 'Chromosomes', chrs.count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Chromosome'})],
            [1, 'Length', chrcontent['length'], 'nt'],
            [1, 'GC-content', ('%0.1f' % (gc_content * 100)), '%'],
        ])
                
        tus = cmodels.TranscriptionUnit.objects.filter(species__id = species.id).annotate(num_genes = Count('genes'))
        mons = tus.filter(num_genes__lte = 1)
        nPolys = tus.filter(num_genes__gt = 1).count()
        content.append([
            [0, 'Transcription units', tus.count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'TranscriptionUnit'})],            
            [1, 'Monocistrons', tus.count() - nPolys],
            [1, 'Polycistrons', nPolys],
        ])
        
        content.append([
            [0, 'Genes', cmodels.Gene.objects.filter(species__id = species.id).count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Gene'})],
            [1, 'mRNA', cmodels.Gene.objects.filter(species__id = species.id, type__wid='mRNA').count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Gene'}) + '?type=mRNA'],
            [1, 'rRNA', cmodels.Gene.objects.filter(species__id = species.id, type__wid='rRNA').count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Gene'}) + '?type=rRNA'],
            [1, 'sRNA', cmodels.Gene.objects.filter(species__id = species.id, type__wid='sRNA').count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Gene'}) + '?type=sRNA'],
            [1, 'tRNA', cmodels.Gene.objects.filter(species__id = species.id, type__wid='tRNA').count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Gene'}) + '?type=tRNA'],
        ])
        
        content.append([
            [0, 'Chromosome features', 
                cmodels.ChromosomeFeature.objects.filter(species__id = species.id).count(),
                None,
                reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'ChromosomeFeature'}),
                ],
            [1, 'DnaA boxes', 
                cmodels.ChromosomeFeature.objects.filter(species__id = species.id, type__parent__wid='ChromosomeFeature-DnaA_box').count(),
                ],                
            [1, 'Short tandem repeats', 
                cmodels.ChromosomeFeature.objects.filter(species__id = species.id, type__parent__wid='ChromosomeFeature-Short_Tandem_Repeat').count(),
                ],
            [1, 'Other', 
                cmodels.ChromosomeFeature.objects.filter(species__id = species.id).exclude(type__parent__wid='ChromosomeFeature-DnaA_box').exclude(type__parent__wid='ChromosomeFeature-Short_Tandem_Repeat').count()],
        ])
        
        content.append([
            [0, 'Metabolites', cmodels.Metabolite.objects.filter(species__id = species.id).count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Metabolite'})],
            [1, 'Amino acids', 
                cmodels.Metabolite.objects.filter(species__id = species.id, type__wid='amino_acid').count() + 
                cmodels.Metabolite.objects.filter(species__id = species.id, type__wid='modified_amino_acid').count() +
                cmodels.Metabolite.objects.filter(species__id = species.id, type__wid='non-standard_amino_acid').count() +
                cmodels.Metabolite.objects.filter(species__id = species.id, type__wid='vitamin_non-standard_amino_acid').count()                
                ],
            [1, 'Antibiotic', 
                cmodels.Metabolite.objects.filter(species__id = species.id, type__wid='antibiotic').count() + 
                cmodels.Metabolite.objects.filter(species__id = species.id, type__parent__wid='antibiotic').count()
                ],
            [1, 'Gases', 
                cmodels.Metabolite.objects.filter(species__id = species.id, type__wid='gas').count() + 
                cmodels.Metabolite.objects.filter(species__id = species.id, type__parent__wid='gas').count()
                ],
            [1, 'Ions', 
                cmodels.Metabolite.objects.filter(species__id = species.id, type__wid='ion').count() + 
                cmodels.Metabolite.objects.filter(species__id = species.id, type__parent__wid='ion').count()
                ],
            [1, 'Lipids', 
                cmodels.Metabolite.objects.filter(species__id = species.id, type__wid='lipid').count() + 
                cmodels.Metabolite.objects.filter(species__id = species.id, type__parent__wid='lipid').count() +
                cmodels.Metabolite.objects.filter(species__id = species.id, type__parent__parent__wid='lipid').count()
                ],
            [1, 'Vitamins', 
                cmodels.Metabolite.objects.filter(species__id = species.id, type__wid='vitamin').count() + 
                cmodels.Metabolite.objects.filter(species__id = species.id, type__parent__wid='vitamin').count()
                ],
        ])
                
        mons = cmodels.ProteinMonomer.objects.filter(species__id = species.id)
        cpxs = cmodels.ProteinComplex.objects.filter(species__id = species.id)
        monDNABind = mons.filter(dna_footprint__length__gt=0).count()
        monIntMem = mons.filter(localization__wid = 'm').exclude(signal_sequence__type = 'lipoprotein').count() + mons.filter(localization__wid = 'tm').exclude(signal_sequence__type = 'lipoprotein').count()
        monLipo = mons.filter(signal_sequence__type = 'lipoprotein').count()
        monSecreted = mons.filter(signal_sequence__type = 'secretory').count()
        monTermOrg = mons.filter(localization__wid = 'tc').count() + mons.filter(localization__wid = 'tm').count()
        cpxDNABind = cpxs.filter(dna_footprint__length__gt=0).count()
        content.append([
            [0, 'Proteins', mons.count() + cpxs.count()],
                [1, 'Monomers', mons.count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'ProteinMonomer'})],            
                    [2, 'DNA-binding', monDNABind],
                    [2, 'Integral membrane', monIntMem],
                    [2, 'Lipoprotein', monLipo, None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'ProteinMonomer'}) + '?signal_sequence__type=lipoprotein'],            
                    [2, 'Secreted', monSecreted, None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'ProteinMonomer'}) + '?signal_sequence__type=secretory'],
                    [2, 'Terminal organelle', monTermOrg],                    
                [1, 'Complexes', cpxs.count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'ProteinComplex'})],
                    [2, 'DNA-binding', cpxDNABind],
        ])

        rxns = cmodels.Reaction.objects.filter(species__id = species.id)
        content.append([
            [0, 'Reactions', rxns.count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'})],
            [1, 'DNA damage', rxns.filter(processes__wid='Process_DNADamage').count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_DNADamage'],
            [1, 'DNA repair', rxns.filter(processes__wid='Process_DNARepair').count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_DNARepair'],
            [1, 'Metabolic', rxns.filter(processes__wid='Process_Metabolism').count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_Metabolism'],            
            [1, 'Protein decay', rxns.filter(processes__wid='Process_ProteinDecay').count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_ProteinDecay'],
            [1, 'Protein modification', rxns.filter(processes__wid='Process_ProteinModification').count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_ProteinModification'],            
            [1, 'Replication Initiation', rxns.filter(processes__wid='Process_ReplicationInitiation').count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_ReplicationInitiation'],
            [1, 'RNA decay', rxns.filter(processes__wid='Process_RNADecay').count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_RNADecay'],
            [1, 'RNA modification', rxns.filter(processes__wid='Process_RNAModification').count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_RNAModification'],            
            [1, 'RNA processing', rxns.filter(processes__wid='Process_RNAProcessing').count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_RNAProcessing'],            
            [1, 'Transcription', rxns.filter(processes__wid='Process_Transcription').count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_Transcription'],            
            [1, 'Translation', rxns.filter(processes__wid='Process_Translation').count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_Translation'],            
            [1, 'tRNA aminoacylation', rxns.filter(processes__wid='Process_tRNAAminoacylation').count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_tRNAAminoacylation'],
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
        
        tr = cmodels.TranscriptionalRegulation.objects.filter(species__id = species.id)
        nTus = len(set([x[0] for x in tr.values_list('transcription_unit')]))
        nTfs = len(set([x[0] for x in tr.values_list('transcription_factor')]))
        content.append([
            [0, 'Transcriptional regulation'],
            [1, 'Interactions', tr.count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'TranscriptionalRegulation'})],
            [1, 'Transcriptional regulators', nTfs],
            [1, 'Regulated promoters', nTus],
        ])        

        content.append([
            [0, 'Pathways', cmodels.Pathway.objects.filter(species__id = species.id).count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Pathway'})],
        ])
        content.append([
            [0, 'Stimuli', cmodels.Stimulus.objects.filter(species__id = species.id).count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Stimulus'})],
        ])
        
        
        nCellComp = cmodels.Metabolite.objects.filter(species__id = species.id, biomass_composition__concentration__isnull=False).count()
        nMediaComp = cmodels.Metabolite.objects.filter(species__id = species.id, media_composition__concentration__isnull=False).count()        
        nKineticsKeq = cmodels.Reaction.objects.filter(species__id = species.id, keq__isnull=False).count()
        nKineticsKm = \
            cmodels.Reaction.objects.filter(species__id = species.id, kinetics_forward__km__isnull=False).count() + \
            cmodels.Reaction.objects.filter(species__id = species.id, kinetics_backward__km__isnull=False).count()
        nKineticsVmax = \
            cmodels.Reaction.objects.filter(species__id = species.id, kinetics_forward__vmax__isnull=False).count() + \
            cmodels.Reaction.objects.filter(species__id = species.id, kinetics_backward__vmax__isnull=False).count()
        nRnaExp = cmodels.Gene.objects.filter(species__id = species.id, expression__isnull=False).count()
        nRnaHl = cmodels.Gene.objects.filter(species__id = species.id, half_life__isnull=False).count()
        nStimuli = cmodels.Stimulus.objects.filter(species__id = species.id, value__isnull=False).count()
        nTrAffinity = cmodels.TranscriptionalRegulation.objects.filter(species__id = species.id, affinity__isnull=False).count()        
        nTrActivity = cmodels.TranscriptionalRegulation.objects.filter(species__id = species.id, activity__isnull=False).count()        
        nOther = cmodels.Parameter.objects.filter(species__id = species.id).count()
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
            [1, 'Other', nOther, None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Parameter'})],            
        ])
        
        content.append([
            [0, 'Processes', cmodels.Process.objects.filter(species__id = species.id).count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Process'})],
        ])
        
        content.append([
            [0, 'States', cmodels.State.objects.filter(species__id = species.id).count(), None, reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'State'})],
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
        #refs = cmodels.PublicationReference.objects.filter(species__id = species.id)
        #sources['total'] = refs.count()
        #sources['types'] = [
        #        {'type': 'Articles', 'count': refs.filter(species__id = species.id, type__wid='article').count()},
        #        {'type': 'Books', 'count': refs.filter(species__id = species.id, type__wid='book').count()},
        #        {'type': 'Thesis', 'count': refs.filter(species__id = species.id, type__wid='thesis').count()},
        #        {'type': 'Other', 'count': refs.filter(species__id = species.id, type__wid='misc').count()},
        #    ]
        #sources['dates'] = refs.filter(year__isnull=False).order_by('year').values('year').annotate(count=Count('year'))
        
            
        nEstimated = cmodels.Parameter.objects.filter(species__id = species.id, value__evidence__is_experimentally_constrained=False).count()
        nExpConstrained = nTotParameters - nEstimated
        sources['evidence_parameters'] = [
            {'type': 'Experimentally constrained', 'count': nExpConstrained},
            {'type': 'Computationally estimated', 'count': nEstimated},
            ]
            
        
        sources['evidence_species'] = cmodels.Evidence.objects.filter(species_component__species__id = species.id).values('species').annotate(count = Count('id'))
        sources['evidence_media'] = cmodels.Evidence.objects.filter(species_component__species__id = species.id).values('media').annotate(count = Count('id'))
        sources['evidence_pH'] = cmodels.Evidence.objects.filter(species_component__species__id = species.id).values('pH').annotate(count = Count('id'))
        sources['evidence_temperature'] = cmodels.Evidence.objects.filter(species_component__species__id = species.id).values('temperature').annotate(count = Count('id'))
            
    return chelpers.render_queryset_to_response(
        species = species,
        data = {
            'content': [contentCol1, contentCol2, contentCol3],
            'contentRows': range(max(len(contentCol1), len(contentCol2), len(contentCol3))),
            'sources': sources,            
            },
        request = request, 
        template = 'cyano/species.html')
    
def about(request, species_wid=None):
    return chelpers.render_queryset_to_response(
        species = species,
        request = request, 
        template = 'cyano/about.html', 
        data = {
            'ROOT_URL': settings.ROOT_URL,
        }
    )    
        
def tutorial(request, species_wid=None):
    return chelpers.render_queryset_to_response(
        species = species,
        request = request, 
        template = 'cyano/tutorial.html')    

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
        species_wid = cmodels.Species.objects.all()[0].wid
    results = SearchQuerySet().filter(species=species).filter(content=query)
    
    #calculate facets        
    facets = results.facet('model_type')
    tmp = facets.facet_counts()['fields']['model_type']
    modelNameFacet = []
    objectTypes = chelpers.getObjectTypes()
    models = []
    for tmp2 in tmp:
        modelName = objectTypes[objectTypes.index(tmp2[0])]
        modelNameFacet.append({
            'name':modelName, 
            'verbose_name': chelpers.getModel(modelName)._meta.verbose_name,
            'count':tmp2[1],
            })
        models.append(chelpers.getModel(modelName))
    modelNameFacet.sort(lambda x, y:cmp(x['verbose_name'], y['verbose_name']))
    
    #narrow search by facets
    model_type = request.GET.get('model_type', '')
    if model_type:
        results = results.models(chelpers.getModel(model_type))
        
    #order results
    results = results.order_by('wid')
    
    #convert results to query set
    queryset = EmptyQuerySet()
    for obj in results:
        tmp = obj.model.objects.none()
        tmp._result_cache.append(obj.object)
        queryset = chain(queryset, tmp)
    
    #form response
    return chelpers.render_queryset_to_response(
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
                
        if isinstance(field, (ForeignKey, ManyToManyField)) and not issubclass(field.rel.to, cmodels.Entry):
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
                id_ = tmp2['wid']
                name = capfirst(tmp2['name'])
            elif (field.choices is not None) and (len(field.choices) > 0) and (not isinstance(field, (BooleanField, NullBooleanField))):    
                id_ = value
                choices = [x[0] for x in field.choices]
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

    if request.is_ajax():
        template = "cyano/list_page.html"
    else:
        template = "cyano/list.html"

    return chelpers.render_queryset_to_response(
        species = species,
        request = request, 
        models = [model], 
        queryset = objects, 
        template = template,
        data = {'facet_fields': facet_fields})

@resolve_to_objects
@permission_required(perm.READ_NORMAL)
def detail(request, species, model, item):
    fieldsets = deepcopy(model._meta.fieldsets)
    
    if request.GET.get('format', 'html') == "html":
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
                    
                data = chelpers.format_field_detail_view(species, item, field_name, request.user.is_anonymous())
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
@permission_required(perm.READ_HISTORY)
def history(request, species, model = None, item = None):
    revisions = []
    entry = []
    date = None

    if item:
        # Item specific
        obj = item
        objects = cmodels.Revision.objects.filter(current = item).distinct("detail").order_by("-detail")
        wid = obj.wid
        detail_id = obj.detail.pk
        date = obj.detail.date.date()
        time = obj.detail.date.strftime("%H:%M")
        reason = obj.detail.reason
        author = obj.detail.user
        url = reverse("cyano.views.detail", kwargs = {"species_wid": species.wid, "model_type": obj.model_type.model_name, "wid": wid})
        
        entry = [date, []]
        entry[1].append({'id': detail_id, 'time': time, 'wid': wid, 'reason': reason, 'author': author, 'url': url})
        
    elif model:
        # Model specific
        tm = cmodels.TableMeta.get_by_model(model)
        components = model.objects.filter(species = species)
        objects = cmodels.Revision.objects.filter(current__model_type = tm, current__pk__in = components).distinct("detail").order_by("-detail")
    else:
        # Whole species specific
        objects = cmodels.Revision.objects.filter(current = item, current__species = species).distinct("detail").order_by("-detail")
    
    for obj in objects:
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
        species = species,
        request = request,
        models = [model],
        queryset = qs,
        template = 'cyano/history.html',
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
    qs = chelpers.objectToQuerySet(item, model = model)

    #render response
    return chelpers.render_queryset_to_response(
        species = species,        
        request = request, 
        models = [model],
        queryset = qs,
        template = 'cyano/history_detail.html', 
        data = {
            'fieldsets': fieldsets,
            'message': request.GET.get('message', ''),
            })
            
@resolve_to_objects
@permission_required(perm.WRITE_NORMAL)
def add(request, species=None, model=None):
    return edit(request, species=species, model=model, action='add')
    
@resolve_to_objects 
@permission_required(perm.WRITE_NORMAL)
def edit(request, species, model = None, item = None, action='edit'):
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
        with transaction.commit_on_success():
            submitted_data = chelpers.get_edit_form_data(model, request.POST, user = request.user.profile)
            
            data = submitted_data
            data['id'] = obj.id
            data['species'] = species.wid
            data['model_type'] = model.__name__
            
            try:
                #validate is WID unique
                if issubclass(model, cmodels.SpeciesComponent):
                    qs = cmodels.SpeciesComponent.objects.values('wid', 'model_type__model_name').filter(species__wid=species.wid)
                else:
                    qs = model.objects.values('wid', 'model_type__model_name').all()
                    
                if action == 'edit':
                    qs = qs.exclude(id=obj.id)
                    
                wids = {}
                for x in qs:
                    wids[x['wid']] = x['model_type__model_name']
                
                if data['wid'] in wids.keys():
                    raise ValidationError({'wid': 'Value must be unique'})
                    
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
                error_messages = error.message_dict

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
            
@resolve_to_objects 
@permission_required(perm.WRITE_DELETE)      
def delete(request, species, model, item):
    qs = chelpers.objectToQuerySet(item, model = model)
    
    #delete
    if request.method == 'POST':
        # Todo: Should be revisioned
        item.delete()
        return HttpResponseRedirect(reverse('cyano.views.listing', kwargs={'species_wid':species.wid, 'model_type': model.__name__}))
        
    #confirmation message
    return chelpers.render_queryset_to_response(
        species = species,
        request = request, 
        models = [model],
        queryset = qs,
        template = 'cyano/delete.html', 
        )

def exportData(request, species_wid=None):
    getDict = request.GET.copy()
    if getDict.get('format', ''):
        getDict.__setitem__('species', getDict.get('species', species_wid))    
    form = ExportDataForm(getDict or None)
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
                queryset = chain(queryset, model.objects.filter(species__id=species.id).select_related(depth=2).all())
            else:
                queryset = chain(queryset, model.objects.select_related(depth=2).filter(id=species.id))
            models.append(chelpers.getModel(model_type))
        
        return chelpers.render_queryset_to_response(
            species = species,
            request = request, 
            queryset = queryset, 
            template = 'cyano/exportDataResult.html', 
            models = models)

@login_required
@resolve_to_objects
def importData(request, species=None):
    if request.method == 'POST':
        form = ImportDataForm(request.POST or None, request.FILES)
        
        if form.is_valid():       
            selected_species_wid = form.cleaned_data['species'] or form.cleaned_data['new_wid']
            if selected_species_wid == '' or selected_species_wid is None:
                form.errors["species"] = ['Please select a species or create a new one']
            elif form.cleaned_data['new_wid'] and cmodels.Species.objects.filter(wid = selected_species_wid).exists():
                form.errors["new_wid"] = ['The identifier specified is already in use']
            else:
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
                request.POST["new_species"] = form.cleaned_data["new_species"]
                request.POST["data_type"] = data_type
                request.POST["reason"] = form.cleaned_data["reason"]
                
                if data_type == "fasta":
                    from cyano.importer.fasta import FastaImporter
                    importer = FastaImporter()
                elif data_type == "fastagene":
                    from cyano.importer.fasta_to_genbank import FastaToGenbankImporter
                    importer = FastaToGenbankImporter()
                elif data_type == "genbank":
                    from cyano.importer.genbank import GenbankImporter
                    importer = GenbankImporter()
                elif data_type == "optgene":
                    from cyano.importer.optgene import OptGeneImporter
                    importer = OptGeneImporter()
                
                if importer:
                    importer.load(filename)
                    return importer.preview(request, selected_species_wid)
                
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

    return chelpers.render_queryset_to_response(
        species = species,
        request = request, 
        template = 'cyano/importDataForm.html', 
        data = {
            'form': form},
        )

@login_required
@resolve_to_objects
def importSubmitData(request, species=None):
    from cyano.importer.genbank import GenbankImporter
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
    errors = chelpers.get_invalid_objects(cmodels.Species.objects.values('id').get(wid=species_wid)['id'])
    
    return chelpers.render_queryset_to_response(
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
    next_url = request.REQUEST.get('next', '')
    
    if request.method == "POST":
        form = AuthenticationForm(data=request.POST)
        if form.is_valid():
            auth_login(request, form.get_user())
            
            if request.session.test_cookie_worked():
                request.session.delete_test_cookie()

            return HttpResponseRedirect(next_url)
    else:
        form = AuthenticationForm(request)

    request.session.set_test_cookie()

    return chelpers.render_queryset_to_response(
        species = species,
        request = request, 
        template = 'cyano/login.html', 
        data = {
            'form': form,
            'next': next_url,
        })

@resolve_to_objects        
def logout(request, species=None):
    auth_logout(request)    
    return chelpers.render_queryset_to_response(
        species = species,
        request = request, 
        template = 'cyano/logout.html', 
        )
    
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
        request = request, 
        template = 'cyano/sitemap_species.xml', 
        data = {
            'ROOT_URL': settings.ROOT_URL,
            'entries': cmodels.SpeciesComponent.objects.filter(species__id = species.id),
        }
    )

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
                    msg = "Invalid POST request on permission page.",
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
            
            # Check if user or group is completly missing but had permissions before
            # Delete if this is the case
            all_user_perms = cmodels.UserPermission.objects.filter(entry = entry).prefetch_related("user")
            for permission in all_user_perms:
                if not permission.user in new_permissions_user_allow and not permission.user in new_permissions_user_deny:
                    permission.delete()
                else:
                    if not permission.user in new_permissions_user_allow:
                        permission.allow.clear()
                    if not permission.user in new_permissions_user_deny:
                        permission.deny.clear()
            
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
                new_perm, _ = cmodels.UserPermission.objects.get_or_create(entry = entry, user = u)
                allows = new_perm.allow.all()
                for i, x in enumerate(range(1, 9)):
                    if p[i] == 0:
                        # delete the entry if it exists
                        if cmodels.Permission.get_by_pk(x) in allows:
                            new_perm.allow.remove(cmodels.Permission.get_by_pk(x))
                    else:
                        # add the entry if missing
                        if not cmodels.Permission.get_by_pk(x) in allows:
                            new_perm.allow.add(cmodels.Permission.get_by_pk(x))
    
            for u, p in new_permissions_user_deny.iteritems():
                new_perm, _ = cmodels.UserPermission.objects.get_or_create(entry = entry, user = u)
                denies = new_perm.deny.all()
                for i, x in enumerate(range(1, 9)):
                    if p[i] == 0:
                        # delete the entry if it exists
                        if cmodels.Permission.get_by_pk(x) in denies:
                            new_perm.deny.remove(cmodels.Permission.get_by_pk(x))
                    else:
                        # add the entry if missing
                        if not cmodels.Permission.get_by_pk(x) in denies:
                            new_perm.deny.add(cmodels.Permission.get_by_pk(x))
    
                new_perm.save()
                        
            for g, p in new_permissions_group_allow.iteritems():
                new_perm, _ = cmodels.GroupPermission.objects.get_or_create(entry = entry, group = g)
                allows = new_perm.allow.all()
                for i, x in enumerate(range(1, 9)):
                    if p[i] == 0:
                        # delete the entry if it exists
                        if cmodels.Permission.get_by_pk(x) in allows:
                            new_perm.allow.remove(cmodels.Permission.get_by_pk(x))
                    else:
                        # add the entry if missing
                        if not cmodels.Permission.get_by_pk(x) in allows:
                            new_perm.allow.add(cmodels.Permission.get_by_pk(x))
                
            for g, p in new_permissions_group_deny.iteritems():     
                new_perm, _ = cmodels.GroupPermission.objects.get_or_create(entry = entry, group = g)    
                denies = new_perm.deny.all()
                for i, x in enumerate(range(1, 9)):
                    if p[i] == 0:
                        # delete the entry if it exists
                        if cmodels.Permission.get_by_pk(x) in denies:
                            new_perm.deny.remove(cmodels.Permission.get_by_pk(x))
                    else:
                        # add the entry if missing
                        if not cmodels.Permission.get_by_pk(x) in denies:
                            new_perm.deny.add(cmodels.Permission.get_by_pk(x))
                
                new_perm.save()

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