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
import json

import math
import os
from crispy_forms.utils import render_crispy_form
from django.contrib.contenttypes.models import ContentType
from django.template.context import Context, RequestContext
from django.views.decorators.csrf import ensure_csrf_cookie
from haystack.inputs import AutoQuery
from jsonview.decorators import json_view
from rest_framework.permissions import IsAuthenticated
import settings
import tempfile
from copy import deepcopy
from itertools import chain

from django.contrib.auth.decorators import login_required, user_passes_test
from django.core.exceptions import ValidationError, ObjectDoesNotExist
from django.core.urlresolvers import reverse
from django.db import transaction
from django.db.models import Count
from django.db.models.fields import BooleanField, NullBooleanField, AutoField, BigIntegerField, DecimalField, FloatField, IntegerField, PositiveIntegerField, PositiveSmallIntegerField, SmallIntegerField
from django.db.models.fields.related import ManyToManyField, ForeignKey
from django.db.models.query import EmptyQuerySet
from django.shortcuts import get_object_or_404
from django.utils.text import capfirst

from haystack.query import SearchQuerySet

from cyano.forms import ExportDataForm, ImportDataForm, ImportSpeciesForm, DeleteForm, CreateBasketForm, \
    RenameBasketForm
import cyano.helpers as chelpers
import cyano.models as cmodels
from cyano.models import PermissionEnum as perm
from cyano.decorators import resolve_to_objects, ajax_required, permission_required
from django.db.transaction import atomic
from django.http.response import HttpResponseRedirect, HttpResponseBadRequest, Http404
from django.http.response import HttpResponse
from rest_framework import generics, filters, permissions
from rest_framework.response import Response
from rest_framework.views import APIView
import cyano.serializers as cserializers


def index(request):
    return chelpers.render_queryset_to_response(
        request,
        template="cyano/index.html",
        )


# Global permissions must be before object resolution
@permission_required(perm.ACCESS_SPECIES)
@resolve_to_objects
@permission_required(perm.READ_NORMAL)
def species(request, species):
    contentcol = []

    for model in chelpers.getModels(cmodels.SpeciesComponent).values():
        contentcol.append(model.get_statistics(species))

    contentcol = filter(lambda x: x is not None, contentcol)

    printcol = []
    for i, x in enumerate(contentcol, start=1):
        for y in x:
            printcol.append([i] + y)

    def chunks(l, n):
        """ Yield successive n-sized chunks from l.
        """
        for i in xrange(0, len(l), n):
            yield l[i:i+n]

    outcol = list(chunks(printcol, int(math.ceil(len(printcol) / 3.0))))

    return chelpers.render_queryset_to_response(
        species=species,
        data={
            'content': outcol,
            'contentRows': range(max(len(p) for p in outcol)),
            },
        request=request,
        template='cyano/species.html')

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

@permission_required(perm.ACCESS_SPECIES)
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

    try:
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
    except TypeError:
        # passing ".." results in TypeError: SWIG director type mismatch in output value of type
        pass

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

@permission_required(perm.ACCESS_SPECIES)
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

    baskets = None
    if request.user.is_authenticated():
        baskets = cmodels.Basket.objects.filter(user=request.user.profile)

    return chelpers.render_queryset_to_response(
        species=species,
        request=request,
        models=[model],
        queryset=objects,
        template=template,
        data={
            'groups': groups,
            'facet_fields': facet_fields,
            'baskets': baskets,
        }
    )


class EntryPermission(permissions.BasePermission):
    """
    Global permission check for blacklisted IPs.
    """

    def has_permission(self, request, view):
        if not request.user.is_authenticated():
            # Special handling for guests
            user = cmodels.UserProfile.objects.get(user__username="guest")
        else:
            user = request.user.profile

        if not user.has_perm(perm.ACCESS_SPECIES):
            return False

        species = view.get_species(view.kwargs["species_wid"])
        return user.has_perms(perm.READ_NORMAL, species)


class EntryList(generics.GenericAPIView):
    filter_backends = (filters.DjangoFilterBackend, filters.SearchFilter,)
    search_fields = ('wid', 'name', 'synonyms')
    filter_fields = ('id', 'wid', 'name')
    filter_class = None
    serializer_class = cserializers.Entry
    permission_classes = (EntryPermission,)

    def get_species(self, species_wid):
        species = cmodels.Species.objects.for_wid(species_wid, get=False)

        try:
            species.get()
        except ObjectDoesNotExist as e:
            raise Http404

        return species

    def get_queryset(self, species_wid, model_type):
        species = self.get_species(species_wid)

        model = chelpers.getModel(model_type)
        if model is None or not issubclass(model, cmodels.SpeciesComponent):
            raise Http404

        return model.objects.for_species(species.get())

    def get(self, request, species_wid, model_type, format=None):
        objects = self.get_queryset(species_wid, model_type)

        model = chelpers.getModel(model_type)
        self.filter_fields += tuple(model._meta.facet_fields)

        self.serializer_class = chelpers.getSerializer(model_type)
        self.filter_class = chelpers.getFilter(model_type)

        objects = self.filter_queryset(objects)

        serializer = chelpers.getSerializer(model_type)(
            objects,
            many=True,
            context={'request': request, 'species': self.get_species(species_wid)})

        data = serializer.data

        return Response(data)

    def post(self, request, species, model, format=None):
        #serializer = EntrySerializer(data=request.DATA)

        #return self.create(request, *args, **kwargs)
        raise Http404


class EntryDetail(APIView):
    serializer_class = cserializers.Entry
    permission_classes = (EntryPermission,)


    def get_species(self, species_wid):
        species = cmodels.Species.objects.for_wid(species_wid, get=False)

        try:
            species.get()
        except ObjectDoesNotExist as e:
            raise Http404

        return species

    def get_object(self, species_wid, model_type, wid):
        species = self.get_species(species_wid)

        try:
            species.get()
        except ObjectDoesNotExist as e:
            raise Http404

        model = chelpers.getModel(model_type)
        if model is None or not issubclass(model, cmodels.SpeciesComponent):
            raise Http404

        species_obj = model.objects.for_species(species.get())

        # resolve numeric wid to real wid
        try:
            i_wid = int(wid)
            obj = species_obj.filter(pk=i_wid)

            if obj.count() != 1:
                raise Http404

            wid = obj[0].pk

        except ValueError:
            pass

        obj = species_obj.for_wid(wid=wid, get=False)

        if obj.count() != 1:
            raise Http404

        return obj

    def get(self, request, species_wid, model_type, wid, format=None):
        objects = self.get_object(species_wid, model_type, wid)
        serializer = chelpers.getSerializer(model_type)(
            objects,
            many=True,
            context={'request': request, 'species': self.get_species(species_wid)}
        )

        return Response(serializer.data)

    def put(self, request, *args, **kwargs):
        return self.update(request, *args, **kwargs)

    def delete(self, request, *args, **kwargs):
        return self.destroy(request, *args, **kwargs)


class BasketList(generics.GenericAPIView):
    filter_backends = (filters.DjangoFilterBackend, filters.SearchFilter,)
    search_fields = ('name',)
    serializer_class = cserializers.Basket
    permission_classes = (IsAuthenticated,)

    def get(self, request):
        objects = cmodels.Basket.objects.filter(user=request.user.profile)
        objects = self.filter_queryset(objects)

        serializer = self.serializer_class(
            objects,
            many=True,
            context={'request': request})

        return Response(serializer.data)

    def post(self, request):
        raise Http404


class BasketDetail(generics.GenericAPIView):
    filter_backends = (filters.DjangoFilterBackend, filters.SearchFilter,)
    serializer_class = cserializers.BasketComponent
    permission_classes = (IsAuthenticated,)

    def get(self, request, basket_id):
        objects = cmodels.BasketComponent.objects.filter(basket__pk=basket_id, basket__user=request.user.profile)
        objects = self.filter_queryset(objects)

        serializer = self.serializer_class(
            objects,
            many=True,
            context={'request': request})

        return Response(serializer.data)

    def post(self, request, basket_id):
        raise Http404


@permission_required(perm.ACCESS_SPECIES)
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
        fieldsets = chelpers.create_detail_fieldset(item, fieldsets, request.user.is_anonymous())

    qs = chelpers.objectToQuerySet(item, model=model)

    baskets = None
    if request.user.is_authenticated():
        baskets = cmodels.Basket.objects.filter(user=request.user.profile)

    #render response
    return chelpers.render_queryset_to_response(
        species=species,
        request=request,
        models=[model],
        queryset=qs,
        template='cyano/detail.html',
        data={
            'fieldsets': fieldsets,
            'message': request.GET.get('message', ''),
            'baskets': baskets
        }
    )


@permission_required(perm.ACCESS_SPECIES)
@resolve_to_objects
@permission_required(perm.READ_NORMAL)
def detail_field(request, species, model, item):
    from django.template import loader
    from django.utils.html import strip_tags

    if request.GET.get('name') is None:
        return HttpResponseBadRequest("Unknown field")

    strip = request.GET.get('strip', False)

    output = chelpers.format_field_detail_view(item, request.GET.get('name'), request.user.is_anonymous())

    if output is None:
        return HttpResponseBadRequest("Unknown field")

    template = loader.get_template("cyano/field.html")
    c = Context({'request': request, 'data': output})

    rendered = template.render(c)

    if strip:
        rendered = strip_tags(rendered)

    return HttpResponse(rendered)


@permission_required(perm.ACCESS_SPECIES)
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


@permission_required(perm.ACCESS_SPECIES)
@resolve_to_objects
@permission_required(perm.READ_HISTORY)
def history_detail(request, species, model, item, detail_id):
    fieldsets = deepcopy(model._meta.fieldsets)
    
    #filter out type, metadata
    fieldsets = filter(lambda x: type(x) == tuple, fieldsets)
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
                if hasattr(field, "get_accessor_name"):
                    verbose_name = field.get_accessor_name()
                else:
                    verbose_name = field.verbose_name
                
            data = chelpers.format_field_detail_view(item, field_name, request.user.is_anonymous(), detail_id)
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


@permission_required(perm.ACCESS_SPECIES)
@resolve_to_objects
@permission_required(perm.WRITE_NORMAL)
def add(request, species=None, model=None):
    return edit(request, species=species, model=model, action='add')


@permission_required(perm.ACCESS_SPECIES)
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
                    return HttpResponseRedirect(obj.get_absolute_url())
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


@permission_required(perm.ACCESS_SPECIES)
@resolve_to_objects 
@permission_required(perm.WRITE_DELETE)
def delete(request, species, model=None, item=None):
    #retrieve object
    if item is None:
        obj = species
    else:
        obj = item
        
    if model is None:
        model = obj.__class__
    
    qs = chelpers.objectToQuerySet(obj, model=model)

    #delete
    if request.method == 'POST':
        # Todo: Should be revisioned with custom message
        rev_detail = cmodels.RevisionDetail(user=request.user.profile, reason="Delete "+obj.wid)
        obj.delete(species, rev_detail)

        if item is None:
            # Delete species
            target_url = reverse('cyano.views.index')
        else:
            target_url = reverse('cyano.views.listing', kwargs={'species_wid': species.wid, 'model_type': model.__name__})

        return HttpResponseRedirect(target_url)

    if item is None:
        # Delete species
        target_url = reverse('cyano.views.delete', kwargs={'species_wid': species.wid})
    else:
        target_url = reverse('cyano.views.delete',
                             kwargs={
                                 'species_wid': species.wid,
                                 'model_type': model.__name__,
                                 'wid': item.wid
                             }
        )

    delete_form = DeleteForm(None)
    delete_form.helper.form_action = target_url

    #confirmation message
    return chelpers.render_queryset_to_response(
        species=species,
        request=request,
        models=[model],
        queryset=qs,
        template='cyano/delete.html',
        data={'delete_form': delete_form}
        )


@permission_required(perm.ACCESS_SPECIES)
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


@permission_required(perm.ACCESS_SPECIES)
@resolve_to_objects
@permission_required(perm.WRITE_NORMAL)
def importData(request, species):
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
        species=species,
        request=request,
        template='cyano/importDataForm.html',
        data=data
    )


@permission_required(perm.ACCESS_SPECIES)
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

            data['success'] = True
            data['message'] = "New %s %s created" % ("mutant" if mutant else "species", species.name)
        else:
            data['success'] = False
            data['message'] = "An error occured. Please check the fields for errors."

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
def login(request, species=None, message=None, error=200, force_next=None):
    from django.contrib.auth.views import login as djlogin

    context = chelpers.get_extra_context(
        species=species,
        request=request,
    )

    context['message'] = message

    if force_next:
        context['next'] = force_next

    response = djlogin(request, extra_context=context)

    if response.status_code != 302:
        response.status_code = error

    return response

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


@permission_required(perm.ACCESS_SPECIES)
@resolve_to_objects
@permission_required(perm.READ_NORMAL)
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

@resolve_to_objects
@permission_required(perm.READ_PERMISSION)
def permission(request, species, model=None, item=None, edit=False):
    from django.contrib.auth.models import User, Group, Permission
    from guardian.shortcuts import get_users_with_perms, get_groups_with_perms, assign_perm, remove_perm
    from itertools import groupby

    if item is not None:
        obj = cmodels.Entry.objects.get(pk=item.pk)
    else:
        obj = cmodels.Entry.objects.get(pk=species.pk)

    perms = map(lambda x: x[0], cmodels.Entry._meta.permissions) + ["change_entry", "delete_entry"]

    permissions = Permission.objects.filter(codename__in=perms)

    def get_permissions():
        up, gp = obj.get_permissions()
        up = up.order_by("user").select_related("user", "permission")
        gp = gp.order_by("group").select_related("group", "permission")

        u = {}
        g = {}

        for k, v in groupby(up, lambda x: x.user):
            u[k] = map(lambda x: x.permission.codename, v)

        for k, v in groupby(gp, lambda x: x.group):
            g[k] = map(lambda x: x.permission.codename, v)

        return u, g

    u, g = get_permissions()

    if request.method == 'POST':
        for r in request.POST:
            if r[:1] == "u":
                try:
                    uid = int(r[1:])
                    usr = User.objects.get(pk=uid)
                    ulist = request.POST[r].strip()
                    if len(ulist) > 0:
                        new_perms = map(lambda x: int(x), ulist.split(" "))
                    else:
                        new_perms = []

                    if usr in u:
                        p = [x for x in permissions.filter(codename__in=u[usr]).values_list("pk", flat=True)]
                    else:
                        p = []

                    for new_perm in new_perms:
                        if new_perm in p:
                            p.remove(new_perm)
                        else:
                            # create new
                            assign_perm(permissions.get(pk=new_perm).codename, usr, obj)

                    for pp in p:
                        # still in list is gone
                        remove_perm(permissions.get(pk=pp).codename, usr, obj)
                except ValueError:
                    continue
            elif r[:1] == "g":
                try:
                    gid = int(r[1:])
                    grp = Group.objects.get(pk=gid)
                    glist = request.POST[r].strip()
                    if len(glist) > 0:
                        new_perms = map(lambda x: int(x), glist.split(" "))
                    else:
                        new_perms = []

                    if grp in g:
                        p = [x for x in permissions.filter(codename__in=g[grp]).values_list("pk", flat=True)]
                    else:
                        p = []

                    for new_perm in new_perms:
                        if new_perm in p:
                            p.remove(new_perm)
                        else:
                            # create new
                            assign_perm(permissions.get(pk=new_perm).codename, grp, obj)

                    for pp in p:
                        # still in list is gone
                        remove_perm(permissions.get(pk=pp).codename, grp, obj)
                except ValueError:
                    continue

    u, g = get_permissions()

    ul = []
    gl = []

    # Add users/groups that aren't
    ux = User.objects.exclude(pk__in=map(lambda x: x.pk, u.keys()))
    gx = Group.objects.exclude(pk__in=map(lambda x: x.pk, g.keys()))

    u.update(dict(zip(ux, [[] for x in range(len(ux))])))
    g.update(dict(zip(gx, [[] for x in range(len(gx))])))

    for usr in u.keys():
        ul.append(
            {"pk": usr.pk,
             "username": usr.username,
             "name": usr.first_name + " " + usr.last_name,
             "permissions": permissions.filter(codename__in=u[usr]).values_list("pk", flat=True)})

    for grp in g.keys():
        gl.append({
            "pk": grp.pk,
            "name": grp.name,
            "permissions": permissions.filter(codename__in=g[grp]).values_list("pk", flat=True)})

    data = {"users": ul, "groups": gl, "permissions": permissions}

    return chelpers.render_queryset_to_response(
        request,
        template="cyano/global_permission.html",
        data=data)

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
    is_admin = request.user.is_superuser

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


@permission_required("access_sbgn")
def sbgn(request):
    return chelpers.render_queryset_to_response(
        request,
        template="cyano/sbgn.html",
    )

@login_required
@ensure_csrf_cookie
@resolve_to_objects
def basket(request, species=None, basket_id=0):
    create_form = CreateBasketForm(None)
    rename_form = RenameBasketForm(None)

    if basket_id == 0:
        bask = None
        template = "cyano/basket.html"

        # Fetch basket list
        queryset = cmodels.Basket.objects.filter(user=request.user.profile).annotate(num_components=Count('components'));
    else:
        template = "cyano/basket_content.html"

        try:
            bask = cmodels.Basket.objects.get(user=request.user.profile, pk=basket_id)
        except ObjectDoesNotExist:
            return chelpers.render_queryset_to_response_error(
                request=request,
                error=404,
                msg="Invalid basket: {}".format(basket_id)
            )

        queryset = bask.components.all().prefetch_related("component")

    return chelpers.render_queryset_to_response(
        species=species,
        request=request,
        queryset=queryset,
        template=template,
        data={'basket': bask, 'create_form': create_form, 'rename_form': rename_form})

@ajax_required
@login_required
@resolve_to_objects
def basket_op(request, species=None):
    from django.http import HttpResponseBadRequest

    # Supported operations:
    # create - Creates a new basket with name 'wid'
    #  return: {id, url}
    # add - Adds item with pk 'wid' to basket 'id'
    #  return {new_op}
    # remove - Removes item with pk 'wid' from basket 'id'
    #  return {new_op}
    # delete - Deletes basket with 'id'
    #  return {}
    # rename - Renames basket with 'id' to 'wid'
    #  return wid
    # list_basket - Lists items in basket 'id'
    #  return {id, items: [list of ids]}
    # list_item - Lists baskets containing item 'id'
    #  return {id, baskets: [list of baskets]}
    # Errors return a bad request

    pk = request.POST.get('id', None)
    wid = request.POST.get('wid', None)
    op = request.POST.get('op', None)
    
    if not op or op not in ["add", "delete", "create", "remove", "rename", "list_basket", "list_item"]:
        return HttpResponseBadRequest("Invalid op")

    if op in ["add", "create", "remove", "rename"] and not wid:
        return HttpResponseBadRequest("Invalid wid")

    if op in ["add", "delete", "remove", "rename"] and not pk:
        return HttpResponseBadRequest("Invalid id")

    if op in ["add", "remove", "list_basket", "list_item"] and not species:
        return HttpResponseBadRequest("Invalid species")

    if op == "create":
        basket = cmodels.Basket.objects.create(user=request.user.profile, name=wid)

        url = {"basket_id": basket.pk}
        if species:
            url.update({"species_wid": species.wid})

        return HttpResponse(json.dumps(
            {"id": basket.pk,
             "url": reverse("cyano-basket", kwargs=url)
            })
        )
    elif op == "delete":
        try:
            cmodels.Basket.objects.get(user=request.user.profile, pk=pk).delete()
            return HttpResponse(json.dumps({}))
        except ObjectDoesNotExist:
            return HttpResponseBadRequest("Invalid basket")
    elif op == "rename":
        try:
            bask = cmodels.Basket.objects.get(user=request.user.profile, pk=pk)
            bask.name = wid
            bask.save()
            return HttpResponse(wid.replace("<", "&lt;").replace(">", "&gt;"))
        except ObjectDoesNotExist:
            return HttpResponseBadRequest("Invalid basket")
    elif op == "list_basket":
        try:
            bask = cmodels.Basket.objects.get(user=request.user.profile, components__species=species.pk, pk=pk)
            items = list(bask.components.filter(species=species.pk).values_list("pk", flat=True))
            return HttpResponse(json.dumps({"id": bask.pk, "items": items}))
        except ObjectDoesNotExist:
            return HttpResponseBadRequest("Invalid basket")
    elif op == "list_item":
        bask = list(cmodels.Basket.objects.filter(user=request.user.profile, components__species=species.pk, components__component=pk).values_list("pk", flat=True))
        return HttpResponse(json.dumps({"id": pk, "baskets": bask}))

    # op is add or remove
    try:
        basket = cmodels.Basket.objects.get(user=request.user.profile, pk=pk)
    except ObjectDoesNotExist:
        return HttpResponseBadRequest("Invalid basket")

    try:
        item = cmodels.SpeciesComponent.objects.for_species(species).get(pk=wid)
    except ObjectDoesNotExist:
        return HttpResponseBadRequest("Invalid item")

    kwargs = {"basket": basket,
              "component": item,
              "species": species}

    if op == "add":
        basket_component, created = cmodels.BasketComponent.objects.get_or_create(**kwargs)
        if not created:
            return HttpResponseBadRequest("Already in basket")

        new_op = "remove"
    else:  # remove
        try:
            cmodels.BasketComponent.objects.get(**kwargs).delete()
        except ObjectDoesNotExist:
            return HttpResponseBadRequest("Not in basket")

        new_op = "add"

    return HttpResponse(json.dumps({"new_op": new_op}))


@login_required
@json_view
def basket_create(request):
    if request.method == 'POST':
        form = CreateBasketForm(request.POST)

        if form.is_valid():
            name = form.cleaned_data.get('name')

            cmodels.Basket.objects.create(user=request.user.profile, name=name)

            return {'success': True}
        else:
            form_html = render_crispy_form(form, context=RequestContext(request))
            return {'success': False, 'form_html': form_html}

    return HttpResponseBadRequest()


@login_required
@json_view
def basket_rename(request, basket_id):
    if request.method == 'POST':
        form = RenameBasketForm(request.POST)

        if form.is_valid():
            name = form.cleaned_data.get('name')

            try:
                basket = cmodels.Basket.objects.get(user=request.user.profile, pk=basket_id)
            except ObjectDoesNotExist:
                return HttpResponseBadRequest()

            basket.name = name
            basket.save()

            return {'success': True}
        else:
            form_html = render_crispy_form(form, context=RequestContext(request))
            return {'success': False, 'form_html': form_html}

    return HttpResponseBadRequest()


@user_passes_test(lambda u: u.is_superuser)
def global_permission(request):
    from django.contrib.auth.models import User, Group, Permission

    u = User.objects.all()
    g = Group.objects.all()

    perms = map(lambda x: x[0], chelpers.get_global_permissions())

    permissions = Permission.objects.filter(codename__in=perms)

    if request.method == 'POST':
        for r in request.POST:
            if r[:1] == "u":
                try:
                    uid = int(r[1:])
                    usr = u.get(pk=uid)
                    usr.user_permissions.remove(*permissions)
                    new_perms = map(lambda x: int(x), request.POST[r].strip().split(" "))
                    usr.user_permissions.add(*new_perms)
                except ValueError:
                    continue
            elif r[:1] == "g":
                try:
                    gid = int(r[1:])
                    grp = g.get(pk=gid)
                    grp.permissions.remove(*permissions)
                    new_perms = map(lambda x: int(x), request.POST[r].strip().split(" "))
                    grp.permissions.add(*new_perms)
                except ValueError:
                    continue

    u = User.objects.prefetch_related("user_permissions")
    g = Group.objects.prefetch_related("permissions")
    ul = []
    gl = []

    for usr in u:
        ul.append(
            {"pk": usr.pk,
             "username": usr.username,
             "name": usr.first_name + " " + usr.last_name,
             "permissions": usr.user_permissions.values_list("pk", flat=True)})

    for grp in g:
        gl.append({
            "pk": grp.pk,
            "name": grp.name,
            "permissions": grp.permissions.values_list("pk", flat=True)
        })

    data = {"users": ul, "groups": gl, "permissions": permissions}

    return chelpers.render_queryset_to_response(
        request,
        template="cyano/global_permission.html",
        data=data)
