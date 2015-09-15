"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""
from django.views.decorators.http import require_POST
from jsonview.decorators import json_view

import kegg.models as models
from kegg.helpers import get_reaction_map, request_extract, items_to_quoted_string
from cyano.decorators import ajax_required, global_permission_required
from cyano.helpers import render_queryset_to_response
from django.core.exceptions import ObjectDoesNotExist
from django.http.response import HttpResponse, HttpResponseBadRequest
from django.contrib.auth.decorators import login_required
from django.db.models.query_utils import Q
from jsonview.exceptions import BadRequest


@global_permission_required("access_kegg")
def index(request):
    items, enzymes, metabolites = request_extract(request)

    ecs = list(map(lambda x: x[0], enzymes))
    metas = list(map(lambda x: x[0], metabolites))

    query = None

    if len(metas) > 0:
        # reduce magic via https://stackoverflow.com/questions/7088173/
        query = reduce(lambda x, y: x | y, [Q(name__icontains=word) for word in metas])
        query |= reduce(lambda x, y: x | y, [Q(title__icontains=word) for word in metas])

    if len(ecs) > 0:
        if query is None:
            query = Q(ec_numbers__name__in=ecs)
        else:
            query |= Q(ec_numbers__name__in=ecs)

    if query is not None:
        kegg_maps = models.Map.objects.filter(query).distinct()
    else:
        kegg_maps = models.Map.objects.all()

    overview = filter(lambda x: x.overview, kegg_maps)
    other = filter(lambda x: not x.overview, kegg_maps)

    return render_queryset_to_response(
        request=request,
        template="kegg/index.html",
        data={
            'overview_maps': overview,
            'other_maps': other,
            'items': items,
            'queries': [] if request.user.is_anonymous() else models.Query.objects.filter(user=request.user.profile),
            'is_anonymous': request.user.is_anonymous()
        }
    )


@global_permission_required("access_kegg")
def map_view(request, map_id):
    export = "export_button" in request.POST

    items, enzymes, metabolites = request_extract(request)
    map_data = get_reaction_map(map_id, enzymes, metabolites, export)

    if export:
        map_data = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>""" + map_data

        import wand.image
        with wand.image.Image(blob=map_data, format="svg") as image:
            png_image = image.make_blob("png")

        response = HttpResponse(png_image, content_type='image/png')
        response['Content-Disposition'] = 'attachment; filename={}.png'.format(map_id)

        return response

    return render_queryset_to_response(
        request=request,
        template="kegg/map.html",
        data={
            'map_id': map_id,
            'map_data': map_data,
            'enzymes': enzymes,
            'metabolites': metabolites,
            'items': items,
            'queries': [] if request.user.is_anonymous() else models.Query.objects.filter(user=request.user.profile),
            'is_anonymous': request.user.is_anonymous()
        }
    )


@require_POST
@ajax_required
@login_required
@global_permission_required("access_kegg")
@json_view
def index_ajax(request):
    pk = request.POST.get("pk", 0)
    op = request.POST.get("op", "")

    if not op in ["load", "delete", "save"]:
        raise BadRequest("Invalid op")

    if op == "save":
        name = request.POST.get("name", "")
        query = request.POST.get("query", "")

        if any(len(x) == 0 for x in [name, query]):
            raise BadRequest("Invalid name or query")

        try:
            query_obj = models.Query.objects.get(name=name, user=request.user.profile)
            query_obj.query = query
            query_obj.save()
            created = False
        except ObjectDoesNotExist:
            query_obj = models.Query.objects.create(name=name, query=query, user=request.user.profile)
            created = True

        return {"name": query_obj.name,
                "pk": query_obj.pk,
                "created": created
        }

    try:
        pk = int(pk)
    except ValueError:
        raise BadRequest("Bad pk")

    try:
        query = models.Query.objects.get(pk=pk, user=request.user.profile)
    except ObjectDoesNotExist:
        raise BadRequest("No item found")

    if op == "load":
        return {"name": query.name, "query": query.query}
    elif op == "delete":
        query.delete()
        return {}

