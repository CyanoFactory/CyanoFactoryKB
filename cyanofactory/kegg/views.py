"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

import kegg.models as models
from kegg.helpers import get_reaction_map, request_extract, items_to_quoted_string
from cyano.decorators import ajax_required, permission_required
from cyano.helpers import render_queryset_to_response
from django.core.exceptions import ObjectDoesNotExist
from django.http.response import HttpResponse, HttpResponseBadRequest
from django.contrib.auth.decorators import login_required
from django.db.models.query_utils import Q


@permission_required("access_kegg")
def index(request):
    items, enzymes, metabolites = request_extract(request)

    ecs = map(lambda x: x[0], enzymes)
    metas = map(lambda x: x[0], metabolites)

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
            'queries': models.Query.objects.filter(user=request.user.profile),
        }
    )


@permission_required("access_kegg")
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
            'queries': models.Query.objects.filter(user=request.user.profile),
        }
    )


@ajax_required
@permission_required("access_kegg")
def index_ajax(request):
    import json

    pk = request.POST.get("pk", 0)
    op = request.POST.get("op", "")

    if not op in ["load", "delete", "save"]:
        return HttpResponseBadRequest("Invalid op")

    if op == "save":
        name = request.POST.get("name", "")
        query = request.POST.get("query", "")

        if all(len(x) == 0 for x in [name, query]):
            return HttpResponseBadRequest("Invalid name or query")

        try:
            query_obj = models.Query.objects.get(name=name, user=request.user.profile)
            query_obj.query = query
            query_obj.save()
            created = False
        except ObjectDoesNotExist:
            query_obj = models.Query.objects.create(name=name, query=query, user=request.user.profile)
            created = True

        return HttpResponse(json.dumps(
            {"name": query_obj.name,
             "pk": query_obj.pk,
             "created": created})
        )

    try:
        pk = int(pk)
    except ValueError:
        return HttpResponseBadRequest("Bad pk")

    try:
        query = models.Query.objects.get(pk=pk, user=request.user.profile)
    except ObjectDoesNotExist:
        return HttpResponseBadRequest("No item found")

    if op == "load":
        return HttpResponse(json.dumps({"name": query.name, "query": query.query}))
    elif op == "delete":
        query.delete()
        return HttpResponse("ok")

