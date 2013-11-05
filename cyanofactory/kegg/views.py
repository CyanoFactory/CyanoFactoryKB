"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""


from django.contrib.auth.decorators import login_required
from django.db.models.query_utils import Q
from cyano.helpers import render_queryset_to_response
import kegg.models as kmodels
from kegg.helpers import get_reaction_map, request_extract, items_to_quoted_string


@login_required
def index(request):
    from itertools import groupby
    from collections import OrderedDict

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
        kegg_maps = kmodels.Map.objects.filter(query).order_by("-overview", "name").distinct()
    else:
        kegg_maps = kmodels.Map.objects.all().order_by("-overview", "name")

    groups = OrderedDict((k, list(v)) for k, v in groupby(kegg_maps, lambda x: x.overview))

    return render_queryset_to_response(
        request=request,
        template="kegg/index.html",
        data={
            'groups': groups,
            'quote': items_to_quoted_string(items),
            'items': items,
        }
    )


@login_required
def map_view(request, map_id):
    items, enzymes, metabolites = request_extract(request)
    map_data = get_reaction_map(map_id, enzymes, metabolites, items)

    return render_queryset_to_response(
        request=request,
        template="kegg/map.html",
        data={
            'map_id': map_id,
            'map_data': map_data,
            'enzymes': enzymes,
            'metabolites': metabolites,
            'items': items,
            'quote': items_to_quoted_string(items),
        }
    )
