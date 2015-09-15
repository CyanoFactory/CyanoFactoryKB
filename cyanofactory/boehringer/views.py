"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Copyright (c) 2013 Roebbe Wuenschiers <wuenschi@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""
from django.template import loader
from django.views.decorators.http import require_POST
from jsonview.decorators import json_view
from jsonview.exceptions import BadRequest

import boehringer.models as models
from boehringer.helpers import format_output
from cyano.decorators import ajax_required, global_permission_required
from cyano.helpers import render_queryset_to_response
from django.core.exceptions import ObjectDoesNotExist
from django.http.response import HttpResponse, HttpResponseBadRequest
from django.shortcuts import render_to_response
from django.template.context import RequestContext
from django.contrib.auth.decorators import login_required


def legacy(request):
    return index(request, True)


@global_permission_required("access_boehringer")
def index(request, legacy=None):
    """
    """
    data = format_output(request)

    if "export_button" in request.POST:
        t = loader.get_template('boehringer/boehringer_svg.html')
        c = RequestContext(request, data)
        response = HttpResponse(t.render(c), content_type='image/svg+xml')
        response['Content-Disposition'] = 'attachment; filename=boehringer.svg'
        return response

    if legacy:
        return render_to_response(
            "boehringer/legacy.html",
            data,
            context_instance=RequestContext(request))
    else:
        return render_queryset_to_response(
            request,
            data=data,
            template="boehringer/index.html"
        )

@require_POST
@ajax_required
@login_required
@global_permission_required("access_boehringer")
@json_view
def index_ajax(request):
    import json

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

        return {
            "name": query_obj.name,
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
