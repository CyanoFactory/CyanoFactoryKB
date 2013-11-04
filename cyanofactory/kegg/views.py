"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""


from django.contrib.auth.decorators import login_required
from cyano.helpers import render_queryset_to_response
import kegg.models as kmodels
from kegg.helpers import get_reaction_map, request_extract_ecs


@login_required
def index(request):
    from itertools import groupby
    from collections import OrderedDict
    from urllib import quote

    enzymes = request_extract_ecs(request)
    items = request.GET.get("items", "")
    items = quote(items)

    ecs = map(lambda x: x[0], enzymes)
    if len(ecs) > 0:
        kegg_maps = kmodels.Map.objects.filter(ec_numbers__name__in=ecs).order_by("-overview", "name")
    else:
        kegg_maps = kmodels.Map.objects.all().order_by("-overview", "name")

    groups = OrderedDict((k, list(v)) for k, v in groupby(kegg_maps, lambda x: x.overview))

    return render_queryset_to_response(
        request=request,
        template="kegg/index.html",
        data={
            'groups': groups,
            'quote': items,
            'items': enzymes,
        }
    )


@login_required
def map_view(request, map_id):
    enzymes = request_extract_ecs(request)
    items = request.GET.get("items", "")
    map_data = get_reaction_map(map_id, enzymes, items)

    return render_queryset_to_response(
        request=request,
        template="kegg/map.html",
        data={
            'map_id': map_id,
            'map_data': map_data,
            'enzymes': enzymes,
            'items': enzymes,
        }
    )

