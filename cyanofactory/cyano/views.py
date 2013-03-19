from copy import deepcopy

# Create your views here.
from django.template import Context, loader
#from polls.models import Choice, Poll
from django.http import HttpResponse, Http404, HttpResponseRedirect, HttpResponseForbidden
from django.core.urlresolvers import reverse
from django.shortcuts import render_to_response, get_object_or_404
from django.template import RequestContext
from django.contrib.auth.forms import AuthenticationForm
from django.db.models.loading import get_model
from django.contrib.auth import login as auth_login, logout as auth_logout

from django.views.decorators.debug import sensitive_post_parameters
from django.views.decorators.cache import never_cache
from django.views.decorators.csrf import csrf_protect

from cyano.helpers import render_queryset_to_response,\
    get_verbose_name_for_field_by_name, render_queryset_to_response_error
from cyano.helpers import get_queryset_object_or_404

from cyano.decorators import read_permission_required

import cyano.models as model

import itertools

def index(request):
    return render_queryset_to_response(
        request,
        template = "cyano/index.html",
        )

@read_permission_required
def species(request, species_wid):
    species = get_queryset_object_or_404(model.Species.objects.filter(wid = species_wid))
    
    #gene_count = len(Seqfeature.filter_by_organism(organism).filter(type_term__name = "gene"))

    contentCol1 = []
    contentCol2 = []
    contentCol3 = []
#===============================================================================

    contentCol1.append([
        [0, 'Genes', model.Gene.objects.filter(species__id = species.id).count, 'nt', reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': 'Gene'})],
    ])
    
    contentCol2.append([
        #[0, 'Chromosomes', 42, None, organism.get_genes_url()],
        [1, 'Length', 100, 'nt'],
        [1, 'GC-content', ('%0.1f' % (0.1 * 100)), '%'],
    ])
    
    contentCol1 = list(itertools.chain.from_iterable(contentCol1))
    contentCol2 = list(itertools.chain.from_iterable(contentCol2))

    return render_queryset_to_response(
                request,
                species = species,
                template = "cyano/species.html",
                data = {
                        "content" : [contentCol1, contentCol2, contentCol3],
                        "contentRows": range(max(len(contentCol1), len(contentCol2), len(contentCol3)))
                       })

def listing(request, species_wid, model_type):
    species = get_queryset_object_or_404(model.Species.objects.filter(wid = species_wid))
    model = get_model("cyano", model_type)
    if model == None:
        raise Http404
    
    objects = model.objects.filter(species__id = species.id)
    
    header = map(lambda x: get_verbose_name_for_field_by_name(model, x), model._meta.listing)
 
    return render_queryset_to_response(
                request,
                species = species,
                template = "cyano/list.html",
                queryset = objects,
                model = model,
                data = {#"contentHeaders" : headers,
                        #"contentRows" : rows,
                        #"content" : content,
                        #"contentRows" : range(len(content[0])),
                        "model_type" : model_type,
                        "header" : header,
                        "field_list" : model._meta.listing,
                })
    pass

def detail(request, species_wid, model_type, wid):
    pass

def permission(request, species_wid = None, model_type = None, wid = None):
    if species_wid != None:
        species = get_queryset_object_or_404(model.Species.objects.filter(wid = species_wid))

    if request.user.is_authenticated():
        # Do something for logged-in users.
        main_group = model.GroupProfile.objects.get(group__name = "Registred")
    else:
        # Do something for anonymous users.
        main_group = model.GroupProfile.objects.get(group__name = "Guest")
    
    
    
    #= request.user.profile
        
        
        return render_queryset_to_response_error(
                request,
                species = species,
                error = 403,
                msg = "You don't have permission to access the permission list for this item")
    
    if model_type != None:
        model = get_model("cyano", model_type)
        if model == None:
            raise Http404
    
    # TODO: Wid
    #entry = get_queryset_object_or_404(species__id = species.id)

    return render_queryset_to_response(
                request,
                species = species,
                #template = Genes.Meta.template,
                data = {"contentHeaders" : contentHeaders,
                        "content" : content,
                        "contentRows" : range(len(content[0])),
                        "urls" : urls,
                })


@sensitive_post_parameters()
@csrf_protect
@never_cache
def login(request, species_wid=None):
    species = None
    if species_wid != None:
        species = get_queryset_object_or_404(model.Species.objects.filter(wid = species_wid))
    
    next_ = request.REQUEST.get('next', '')
    
    if request.method == "POST":
        form = AuthenticationForm(data=request.POST)
        if form.is_valid():
            auth_login(request, form.get_user())
            
            if request.session.test_cookie_worked():
                request.session.delete_test_cookie()

            return HttpResponseRedirect(next_)
    else:
        form = AuthenticationForm(request)

    request.session.set_test_cookie()

    return render_queryset_to_response(
        species = species,
        request = request, 
        template = 'cyano/login.html', 
        data = {
            'form': form,
            'next': next_,
        })
        
def logout(request, species_wid=None):
    species = None
    if species_wid != None:
        species = get_queryset_object_or_404(model.Species.objects.filter(wid = species_wid))

    auth_logout(request)    
    return render_queryset_to_response(
        species = species,
        request = request, 
        template = 'cyano/logout.html', 
        )

def genes(request, organism):
    organism = get_queryset_object_or_404(Bioentry.objects.filter(biodatabase__name = database_name).filter(name = organism))
    record = get_database_item(organism)

    contentHeaders = [x[1] for x in Genes.Meta.field_list]
    content = [[] for x in contentHeaders]
    urls = []
    
    genes = []
    cds = {}
    for item in record.features:
        if item.type == "gene":
            if "gene" in item.qualifiers:
                genes.append(item.qualifiers["gene"][0])
        if item.type == "CDS":
            if "gene" in item.qualifiers:
                cds[item.qualifiers["gene"][0]] = item
    
    for item in genes:
        #print item
        if item in cds:
            val = cds[item].qualifiers
            i = 0
            #for x in Genes.Meta.field_list:
            #    content[i].append(val.get(x[0], [""])[0])
            #    i += 1
            #    if (x[0] == Genes.Meta.field_url):
            #        urls.append(organism.get_gene_url(val.get(x[0], [""])[0]))
    
    #for item in itertools.ifilter(lambda x : x.type == "gene", record):
    #    if "gene" in item.qualifiers:
    #        contentCol1.append(item.qualifiers["gene"][0])
    #        contentCol2.append(item.qualifiers["locus_tag"][0])
    #        urls.append(organism.get_gene_url(contentCol1[-1]))
    #    else:
    #        print "Missing"
    #        print item.qualifiers

    return render_queryset_to_response(
                request,
                organism = organism,
                #template = Genes.Meta.template,
                data = {"contentHeaders" : contentHeaders,
                        "content" : content,
                        "contentRows" : range(len(content[0])),
                        "urls" : urls,
                })
    
def gene(request, organism, gene):
    organism = get_queryset_object_or_404(Bioentry.objects.filter(biodatabase__name = database_name).filter(name = organism))
    gene_name = gene
    gene = Gene(organism, gene)
    
    
    
    print gene.get_databases_html()

    return render_queryset_to_response(
                request,
                organism = organism,
                #template = Gene.Meta.template,
                data = {"name": gene_name,
                        "type": "Gene",
                        "fieldsets": gene.fieldsets}
                )
