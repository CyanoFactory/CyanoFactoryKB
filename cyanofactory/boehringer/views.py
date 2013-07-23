from django.shortcuts import render_to_response
import boehringer.models as models
from django.template.context import RequestContext
import re
from cyano.helpers import render_queryset_to_response
from django.contrib.auth.decorators import login_required

def legacy(request):
    return index(request, True)

@login_required
def index(request, legacy = None):
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
