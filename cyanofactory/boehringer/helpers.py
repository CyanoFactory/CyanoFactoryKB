import base64
import re
from django.conf import settings
from django.contrib.staticfiles import finders
import boehringer.models as models

def format_output(request):
    if not "items" in request.POST:
        items = [["1.1.1.1", None, None],
                 ["2.2.2.2", None, None],
                 ["4.1.2.20", "hsl(240.0,100%,50%)", 0],
                 ["1.1.1.20", "hsl(0.0,100%,50%)", 100],
                 ["1.2.1.3", "green", "green"],
                 ["ascorbate", "yellow", "yellow"]]
    else:
        items = []
        for item in re.split("\s+", request.POST["items"]):
            if len(item) == 0:
                continue

            if "#" in item:
                item = item.split("#", 2)
                if len(item[0]) == 0:
                    continue
                item.append(item[1])
                try:
                    h = (100 - int(item[2])) * 2.4
                    item[2] = "hsl(" + str(h) + ", 100%, 50%)"
                except ValueError:
                    pass
                items.append([item[0], item[2], item[1]])
            else:
                items.append([item, None, None])
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
        hits = models.Metabolite.objects.filter(title__icontains=metabolite[0]).prefetch_related("color")
        if hits.count() > 0:
            metabolites_hits += 1
            for x in hits:
                metabolites.append([x, metabolite[1]])
    enzymes_hits = 0
    for enzyme in enzyme_items:
        hits = models.Enzyme.objects.filter(ec=enzyme[0]).prefetch_related("color")
        if hits.count() > 0:
            enzymes_hits += 1
            for x in hits:
                enzymes.append([x, enzyme[1]])

    finder = finders.find("boehringer/pathways.png")

    export = "export_button" in request.POST

    if export:
        with open(finder, "rb") as image_file:
            encoded_string = base64.b64encode(image_file.read())
    else:
        encoded_string = ""

    data = {'image': "data:image/png;base64," + encoded_string.decode("utf-8"),
            'export': export,
            'items': items,
            'metabolites': metabolites,
            'enzymes': enzymes,
            'metabolites_hits': metabolites_hits,
            'enzymes_hits': enzymes_hits,
            'metabolites_no_hits': len(metabolite_items) - metabolites_hits,
            'enzymes_no_hits': len(enzyme_items) - enzymes_hits,
            'queries': [] if request.user.is_anonymous() else models.Query.objects.filter(user=request.user.profile),
            'is_anonymous': request.user.is_anonymous()
    }

    return data
