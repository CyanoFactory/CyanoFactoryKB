import json
from django.http.response import HttpResponse
from cyano.decorators import ajax_required
from PyNetMet.metabolism import *
from PyNetMet.enzyme import *
from PyNetMet.fba import *
from django.conf import settings
from cyano.helpers import render_queryset_to_response


def index(request):
    data = {}
    
    return render_queryset_to_response(request, template="cyanodesign/index.html", data=data)


@ajax_required
def get_reactions(request):
    model = "{}/cyanodesign/models/{}".format(settings.ROOT_DIR, "toy_model.txt")

    org = Metabolism(model)

    ret = []

    for enzyme, constr in zip(org.enzymes, org.constr):
        ret.append({
            "name": enzyme.name,
            "stoichiometric": enzyme.stoic,
            "substrates": enzyme.substrates,
            "products": enzyme.products,
            "reversible": enzyme.reversible,
            "constraints": constr
        })

    return HttpResponse(json.dumps({"external": org.external, "enzymes": ret}), content_type="application/json")


@ajax_required
def calculate(request):
    enzymes = []



    # POST: reactions, constraints, external

    for reaction in reactions:
        enzymes.append(Enzyme(reaction))


@ajax_required
def simulate(request):
    data = json.loads(request.GET["enzymes"])
    constr = json.loads(request.GET["constraints"])
    ext = json.loads(request.GET["external"])
    objective = json.loads(request.GET["objective"])

    print data
    print constr
    print ext
    print objective

    organism = Metabolism("model_name",
                          reactions=data,
                          constraints=constr,
                          external=ext,
                          objective=objective,
                          fromfile=False)

    fba = FBA(organism)

    return HttpResponse(fba)


def export(request):
    pass
