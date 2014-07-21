"""
Copyright (c) 2014 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

import StringIO
import json
import os
import tempfile
from django.contrib.auth.decorators import login_required
from django.core.exceptions import ObjectDoesNotExist
from django.http.response import HttpResponse, HttpResponseBadRequest
from django.shortcuts import redirect
from django.views.decorators.csrf import ensure_csrf_cookie
from cyano.decorators import ajax_required
from PyNetMet.metabolism import *
from PyNetMet.fba import *
from cyano.helpers import render_queryset_to_response, render_queryset_to_response_error
from cyanodesign.forms import UploadModelForm
from cyanodesign.helpers import model_from_string, apply_commandlist
from .models import DesignModel

@login_required
def index(request):
    models = DesignModel.objects.filter(user=request.user.profile)

    return render_queryset_to_response(
        request,
        template="cyanodesign/list.html",
        queryset=models,
        data={}
    )


@login_required
@ensure_csrf_cookie
def design(request, pk):
    try:
        item = DesignModel.objects.get(user=request.user.profile, pk=pk)
    except ObjectDoesNotExist:
        return render_queryset_to_response_error(request, error=404, msg="Model not found")

    data = {"pk": pk, "name": item.name}
    
    return render_queryset_to_response(request, template="cyanodesign/design.html", data=data)


@login_required
@ajax_required
def get_reactions(request, pk):
    #model = "{}/cyanodesign/models/{}".format(settings.ROOT_DIR, "toy_model.txt")

    try:
        item = DesignModel.objects.get(user=request.user.profile, pk=pk).content
    except:
        return HttpResponseBadRequest("Bad Model")

    org = model_from_string(item)

    ret = []

    for enzyme, constr in zip(org.enzymes, org.constr):
        # Filter out auto-created external transport funcs
        if not enzyme.name.endswith("_ext_transp"):
            ret.append({
                "name": enzyme.name,
                "stoichiometric": enzyme.stoic,
                "substrates": enzyme.substrates,
                "products": enzyme.products,
                "reversible": enzyme.reversible,
                "constraints": constr
            })

    return HttpResponse(json.dumps({"external": org.external, "enzymes": ret, "objective": org.objective}), content_type="application/json")


@login_required
@ajax_required
def simulate(request, pk):
    if not all(x in request.POST for x in ["commands", "disabled", "objective"]):
        return HttpResponseBadRequest("Request incomplete")

    try:
        content = DesignModel.objects.get(user=request.user.profile, pk=pk).content
    except ObjectDoesNotExist:
        return HttpResponseBadRequest("Bad Model")

    org = model_from_string(content)

    try:
        commands = json.loads(request.POST["commands"])
        disabled = json.loads(request.POST["disabled"])
        objective = json.loads(request.POST["objective"])
    except ValueError:
        return HttpResponseBadRequest("Invalid JSON data")

    if not all(isinstance(x, list) for x in [commands, disabled, objective]):
        return HttpResponseBadRequest("Invalid data type")

    string_type = map(lambda x: all(isinstance(y, basestring) for y in x), [commands, disabled, objective])
    if not all(x for x in string_type):
        return HttpResponseBadRequest("Invalid data type")


    org = apply_commandlist(org, commands)

    #print "\n".join(data)
    #print constr
    #print ext
    #print objective

    try:
        fba = FBA(org)
    except ValueError:
        return HttpResponseBadRequest("FBA error")

    return HttpResponse(fba)


@login_required
@ajax_required
def export(request):
    pass

@login_required
def upload(request):
    data = {}

    if request.method == 'POST':
        form = UploadModelForm(request.POST, request.FILES)

        if form.is_valid():
            name = form.cleaned_data.get('name')

            #save to temporary file
            filename = request.FILES['file'].name

            ss = StringIO.StringIO()

            with tempfile.NamedTemporaryFile(delete=False) as fid:
                path = fid.name

                for chunk in request.FILES['file'].chunks():
                    ss.write(chunk)
                    fid.write(chunk)

            try:
                Metabolism(path)
                os.remove(path)
            except:
                os.remove(path)
                return HttpResponseBadRequest("Bad Model")

            DesignModel.objects.create(
                user=request.user.profile,
                name=name,
                filename=filename,
                content=ss.getvalue()
            )

            return redirect("cyano-design-index")

    return HttpResponseBadRequest()

@login_required
@ajax_required
def delete(request):
    pk = request.POST.get("id", 0)

    try:
        DesignModel.objects.get(user=request.user.profile, pk=pk).delete()
    except:
        return HttpResponseBadRequest("Bad Model")

    return HttpResponse("ok")
