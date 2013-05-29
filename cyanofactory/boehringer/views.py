from django.http import HttpResponseBadRequest
from django.shortcuts import render_to_response
from django.shortcuts import redirect
from django.template import Context, loader
from db_xref.helpers import get_database_url
import boehringer.models as models
from django.template.context import RequestContext
from django.core.exceptions import ObjectDoesNotExist
import re

def index(request):
    """
    """
    #models.BioMolecule
    if not "items" in request.GET:
        items = [["1.1.1.1", None],
                 ["2.2.2.2", None],
                 ["4.1.2.20", "green"],
                 ["1.1.1.20", None],
                 ["1.2.1.3", None],
                 ["ascorbate", "red"]]
    else:
        items = []
        for item in re.split("\s+", request.GET["items"]):
            if len(item) == 0:
                continue

            if "#" in item:
                item = item.split("#", 2)
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
        hits = models.Metabolite.objects.filter(title__contains = metabolite[0]).prefetch_related("color")
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
    
    print "enzymes"
    for x in enzyme_items:
        print x
    
    print "metabolites"
    for x in metabolite_items:
        print x
    
    data = {}
    data['items'] = items
    data['metabolites'] = metabolites
    data['enzymes'] = enzymes
    data['metabolites_hits'] = metabolites_hits
    data['enzymes_hits'] = enzymes_hits
    data['metabolites_no_hits'] = len(metabolite_items) - metabolites_hits
    data['enzymes_no_hits'] = len(enzyme_items) - enzymes_hits

    return render_to_response(
        "boehringer/index.html",
        data,
        context_instance = RequestContext(request))
