"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Copyright (c) 2013 Roebbe Wuenschiers <wuenschi@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

import re
import boehringer.models as models
from cyano.decorators import ajax_required
from cyano.helpers import render_queryset_to_response
from django.core.exceptions import ObjectDoesNotExist
from django.http.response import HttpResponse, HttpResponseBadRequest
from django.shortcuts import render_to_response
from django.template.context import RequestContext
from django.contrib.auth.decorators import login_required


def legacy(request):
    return index(request, True)


@login_required
def index(request, legacy=None):
    """
    """
    if not "items" in request.POST:
        items = [["1.1.1.1", None],
                 ["2.2.2.2", None],
                 ["4.1.2.20", "green"],
                 ["1.1.1.20", None],
                 ["1.2.1.3", None],
                 ["ascorbate", "red"]]
    else:
        items = []
        for item in re.split("\s+", request.POST["items"]):
            if len(item) == 0:
                continue

            if "#" in item:
                item = item.split("#", 2)
                if len(item[0]) == 0:
                    continue
                items.append([item[0], item[1]])
            else:
                items.append([item, None])
    
    metabolite_items = []
    enzyme_items = []
    
    for item in items:
        maybe_enzyme = item[0].split(".")
        if len(maybe_enzyme) == 4:
            try:
                map(lambda x: int(x), maybe_enzyme)
                enzyme_items.append(item)
                continue
            except ValueError:
                # not a valid EC number, maybe a metabolite
                pass

        metabolite_items.append(item)
    
    metabolites = []
    enzymes = []
    metabolites_hits = 0
    for metabolite in metabolite_items:
        hits = models.Metabolite.objects.filter(title__icontains = metabolite[0]).prefetch_related("color")
        if hits.count() > 0:
            metabolites_hits += 1
            for x in hits:
                metabolites.append([x, metabolite[1]])
    
    enzymes_hits = 0
    for enzyme in enzyme_items:
        hits = models.Enzyme.objects.filter(ec = enzyme[0]).prefetch_related("color")
        if hits.count() > 0:
            enzymes_hits += 1
            for x in hits:
                enzymes.append([x, enzyme[1]])

    data = {}
    data['items'] = items
    data['metabolites'] = metabolites
    data['enzymes'] = enzymes
    data['metabolites_hits'] = metabolites_hits
    data['enzymes_hits'] = enzymes_hits
    data['metabolites_no_hits'] = len(metabolite_items) - metabolites_hits
    data['enzymes_no_hits'] = len(enzyme_items) - enzymes_hits
    data['queries'] = models.Query.objects.filter(user=request.user.profile)

    if legacy:
        return render_to_response(
            "boehringer/legacy.html",
            data,
            context_instance = RequestContext(request))
    else:
        return render_queryset_to_response(
            request,
            data = data,
            template = "boehringer/index.html"
        )


@ajax_required
@login_required
def index_ajax(request):
    import json

    pk = request.GET.get("pk", 0)
    op = request.GET.get("op", "")

    if not op in ["load", "delete", "save"]:
        return HttpResponseBadRequest("Invalid op")

    if op == "save":
        name = request.GET.get("name", "")
        query = request.GET.get("query", "")

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

