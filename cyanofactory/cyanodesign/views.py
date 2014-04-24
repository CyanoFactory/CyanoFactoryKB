import json
from django.contrib.auth.decorators import login_required
from django.http.response import HttpResponse, HttpResponseBadRequest
from cyano.decorators import ajax_required
from PyNetMet.metabolism import *
from PyNetMet.fba import *
from django.conf import settings
from cyano.helpers import render_queryset_to_response


@login_required
def index(request):
    data = {}
    
    return render_queryset_to_response(request, template="cyanodesign/index.html", data=data)


@login_required
@ajax_required
def get_reactions(request):
    model = "{}/cyanodesign/models/{}".format(settings.ROOT_DIR, "toy_model2.txt")

    org = Metabolism(model)

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
def simulate(request):
    if not all(x in request.GET for x in ["enzymes", "constraints", "external", "objective"]):
        return HttpResponseBadRequest("Request incomplete")

    try:
        data = json.loads(request.GET["enzymes"])
        constr = json.loads(request.GET["constraints"])
        ext = json.loads(request.GET["external"])
        objective = json.loads(request.GET["objective"])
    except ValueError:
        return HttpResponseBadRequest("Invalid JSON data")

    if not all(isinstance(x, list) for x in [data, constr, ext, objective]):
        return HttpResponseBadRequest("Invalid data type")

    string_type = map(lambda x: all(isinstance(y, basestring) for y in x), [data, constr, ext, objective])
    if not all(x for x in string_type):
        return HttpResponseBadRequest("Invalid data type")

    #print "\n".join(data)
    #print constr
    #print ext
    #print objective

    try:
        organism = Metabolism("model_name",
                              reactions=data,
                              constraints=constr,
                              external=ext,
                              objective=objective,
                              fromfile=False)

        fba = FBA(organism)
    except ValueError:
        return HttpResponseBadRequest("Invalid data type")

    return HttpResponse(fba)


def export(request):
    pass
