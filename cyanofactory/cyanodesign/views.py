"""
Copyright (c) 2014 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

import StringIO
from collections import OrderedDict
import os
import tempfile
from crispy_forms.utils import render_crispy_form
from django.core.exceptions import ObjectDoesNotExist
from django.http.response import HttpResponse, HttpResponseBadRequest
from django.template.context import RequestContext
from django.views.decorators.csrf import ensure_csrf_cookie
from jsonview.decorators import json_view
from networkx.algorithms.shortest_paths.generic import has_path, all_shortest_paths
from networkx.exception import NetworkXNoPath
from cyano.decorators import ajax_required, permission_required
from PyNetMet2.metabolism import Metabolism
from cyano.helpers import render_queryset_to_response, render_queryset_to_response_error
from cyanodesign.forms import UploadModelForm
from cyanodesign.helpers import model_from_string, apply_commandlist
from .models import DesignModel
import networkx as nx
from networkx.readwrite import json_graph
import math
import pygraphviz
import json


@permission_required("access_cyanodesign")
def index(request):
    upload_form = UploadModelForm(None)

    models = DesignModel.objects.filter(user=request.user.profile)

    return render_queryset_to_response(
        request,
        template="cyanodesign/list.html",
        queryset=models,
        data={'upload_form': upload_form}
    )


@ensure_csrf_cookie
@permission_required("access_cyanodesign")
def design(request, pk):
    try:
        item = DesignModel.objects.get(user=request.user.profile, pk=pk)
    except ObjectDoesNotExist:
        return render_queryset_to_response_error(request, error=404, msg="Model not found")

    data = {"pk": pk, "name": item.name}

    return render_queryset_to_response(request, template="cyanodesign/design.html", data=data)


@ajax_required
@permission_required("access_cyanodesign")
def get_reactions(request, pk):
    try:
        item = DesignModel.objects.get(user=request.user.profile, pk=pk).content
    except:
        return HttpResponseBadRequest("Bad Model")

    org = model_from_string(item)

    if org is None:
        return HttpResponseBadRequest("Bad Model")

    ret = []

    for enzyme in org.enzymes:
        # Filter out auto-created external transport funcs
        if not enzyme.name.endswith("_transp"):
            ret.append({
                "name": enzyme.name,
                "stoichiometric": enzyme.stoic,
                "substrates": enzyme.substrates,
                "products": enzyme.products,
                "reversible": enzyme.reversible,
                "constraints": enzyme.constraint,
                "disabled": enzyme.disabled
            })

    outgraph = ""

    return HttpResponse(json.dumps({
        "external": filter(lambda x: not x.startswith("#"), org.external),
        "enzymes": ret,
        "objective": org.obj[0] if org.obj else None,
        "graph": outgraph}), content_type="application/json")


@permission_required("access_cyanodesign")
def simulate(request, pk):
    if not all(x in request.POST for x in ["commands", "disabled", "objective", "display", "auto_flux"]):
        return HttpResponseBadRequest("Request incomplete")

    output_format = request.POST.get("format", "json")
    dry_run = request.POST.get("dry_run", False)

    try:
        model = DesignModel.objects.get(user=request.user.profile, pk=pk)
        content = model.content
    except ObjectDoesNotExist:
        return HttpResponseBadRequest("Bad Model")

    org = model_from_string(content)

    try:
        commands = json.loads(request.POST["commands"])
        disabled = json.loads(request.POST["disabled"])
        objective = json.loads(request.POST["objective"])
        display = json.loads(request.POST["display"])
        auto_flux = json.loads(request.POST["auto_flux"])
    except ValueError:
        return HttpResponseBadRequest("Invalid JSON data")

    try:
        auto_flux = bool(auto_flux)
    except ValueError:
        return HttpResponseBadRequest("Invalid data type")

    if not all(isinstance(x, list) for x in [commands, disabled, display]):
        return HttpResponseBadRequest("Invalid data type")

    if not isinstance(objective, basestring):
        return HttpResponseBadRequest("Invalid data type")

    try:
        org = apply_commandlist(org, commands)
        if not objective:
            raise ValueError("No objective specified")
        obj_reac = org.get_reaction(objective)
        if obj_reac is None:
            raise ValueError("Objective not in model: " + objective)
        org.obj = [[obj_reac.name, "1"]]
    except ValueError as e:
        return HttpResponseBadRequest("Model error: " + e.message)

    for reaction in org.enzymes:
        reaction.disabled = False

    for item in disabled:
        try:
            reac = org.get_reaction(item)
        except ValueError:
            return HttpResponseBadRequest("Bad disabled list: " + item)

        reac.disabled = True

    if obj_reac.disabled:
        return HttpResponseBadRequest("Objective disabled: " + objective)

    display = filter(lambda x: not len(x) == 0, display)

    for item in display:
        if not org.has_reaction(item):
            return HttpResponseBadRequest("Bad display list: " + item)

    try:
        fba = org.fba()
    except ValueError as e:
        return HttpResponseBadRequest("FBA error: " + e.message)

    if dry_run:
        return HttpResponse(
            json.dumps({"solution": org.solution}),
            content_type="application/json")

    full_g, nodeIDs = calcReactions(org, fba)

    import itertools

    display = json.loads(request.POST["display"])

    # Auto filter by FBA
    if auto_flux:
        flux = filter(lambda x: not x[0].endswith("_transp"), zip(map(lambda x: x.name, fba.reacs), fba.flux))
        flux = sorted(flux, key=lambda x: -x[1])
        display = map(lambda x: x[0], flux[:20])

    g = getSelectedReaction(json_graph.node_link_data(full_g), nodeIDs, display, org)

    graph = nx.to_agraph(g)
    graph.graph_attr.update(splines=True, overlap=False, rankdir="LR")
    graph.node_attr.update(style="filled", colorscheme="pastel19")
    outgraph = str(graph.to_string())
    outgraph = pygraphviz.AGraph(outgraph)
    outgraph = outgraph.draw(format="svg", prog="dot")

    if output_format == "json":
        return HttpResponse(json.dumps(
            {"graph": outgraph,
            "solution": fba.get_status()
            }),
            content_type="application/json"
        )
    elif output_format == "png":
        import wand.image
        with wand.image.Image(blob=outgraph, format="svg") as image:
            png_image = image.make_blob("png")

        r = HttpResponse(png_image, content_type="image/png")
        r['Content-Disposition'] = 'attachment; filename={}.png'.format(model.filename)
        return r
    elif output_format == "svg":
        r = HttpResponse(outgraph, content_type="image/svg+xml")
        r['Content-Disposition'] = 'attachment; filename={}.svg'.format(model.filename)
        return r
    elif output_format == "csv":
        s = StringIO.StringIO()
        for reac in org.enzymes:
            s.write(reac.name)
            s.write("\t")
            s.write(reac.flux)
            s.write("\r\n")
        r = HttpResponse(s.getvalue(), content_type="text/csv")
        r['Content-Disposition'] = 'attachment; filename={}.csv'.format(model.filename)
        return r
    else:
        return HttpResponseBadRequest("Unknown format")


@permission_required("access_cyanodesign")
def export(request, pk):
    try:
        model = DesignModel.objects.get(user=request.user.profile, pk=pk)
        content = model.content
    except ObjectDoesNotExist:
        return HttpResponseBadRequest("Bad Model")

    response = HttpResponse(
        content,
        mimetype="application/x-bioopt; charset=UTF-8",
        content_type="application/x-bioopt; charset=UTF-8"
    )
    response['Content-Disposition'] = "attachment; filename=" + model.filename

    return response


@ajax_required
@permission_required("access_cyanodesign")
def save(request, pk):
    if not all(x in request.POST for x in ["commands", "disabled", "objective"]):
        return HttpResponseBadRequest("Request incomplete")

    try:
        model = DesignModel.objects.get(user=request.user.profile, pk=pk)
        content = model.content
    except ObjectDoesNotExist:
        return HttpResponseBadRequest("Bad Model")

    org = model_from_string(content)

    try:
        commands = json.loads(request.POST["commands"])
        disabled = json.loads(request.POST["disabled"])
        objective = json.loads(request.POST["objective"])
    except ValueError:
        return HttpResponseBadRequest("Invalid JSON data")

    if not all(isinstance(x, list) for x in [commands, disabled]):
        return HttpResponseBadRequest("Invalid data type")

    if not isinstance(objective, basestring):
        return HttpResponseBadRequest("Invalid data type")

    try:
        org = apply_commandlist(org, commands)
        if objective:
            obj_reac = org.get_reaction(objective)
            if obj_reac is None:
                raise ValueError("Objective not in model: " + objective)
            org.objective = obj_reac
        else:
            org.objective = None
    except ValueError as e:
        return HttpResponseBadRequest("Model error: " + e.message)

    for reaction in org.enzymes:
        reaction.disabled = False

    for item in disabled:
        try:
            reac = org.get_reaction(item)
        except ValueError:
            return HttpResponseBadRequest("Bad disabled list: " + item)

        reac.disabled = True

    from StringIO import StringIO
    sio = StringIO()
    org.write(sio)
    model.content = sio.getvalue()
    model.save()

    return HttpResponse("OK")


@permission_required("access_cyanodesign")
@json_view
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

            return {'success': True}
        else:
            form_html = render_crispy_form(form, context=RequestContext(request))
            return {'success': False, 'form_html': form_html}

    return HttpResponseBadRequest()


@permission_required("access_cyanodesign")
@ajax_required
def delete(request):
    pk = request.POST.get("id", 0)

    try:
        DesignModel.objects.get(user=request.user.profile, pk=pk).delete()
    except:
        return HttpResponseBadRequest("Bad Model")

    return HttpResponse("ok")

def flatten(li):
    return [item for sublist in li for item in sublist]

def getSelectedReaction(jsonGraph, nodeDic, reacIDs, org):
    """
    Filtering selected Reactions and show Results from PyNetMet calculation.
    It returns a subgraph of the Graph from the jsonGraph. The output is a DOT-Language String.
    @param jsonGraph: Graph in JSON-Format
    @param nodeDic: dict mapping names to ids
    @param reacIDs: Name of reactions that contained in the nodeDic
    @param org: organism
    @return Subgraph
    """
    # Translate reac names to IDs
    # Get substrates and products of all reacs
    metabolites = []
    for reac in reacIDs:
        metabolites += org.get_reaction(reac).metabolites

    met_ids = map(lambda x: nodeDic[x], metabolites)
    g = json_graph.node_link_graph(jsonGraph)

    g.remove_edges_from(filter(lambda x: g.get_edge_data(*x)["object"].name not in reacIDs, g.edges(met_ids)))

    # Get products/substrates directly connected to filter
    #reacIDs += flatten(g.in_edges(reacIDs)) + flatten(g.out_edges(reacIDs))

    h = g.subgraph(met_ids)
    return h

def calcReactions(org, fluxResults):
    """
    Adding Results from FLUX-Analysis
    @param org: organism
    @param fluxResults: Results from PyNetMet FLUX-Analysis Calculation as String
    @return Graph with added Attributes
    """
    changeDic = {}
    reactions = fluxResults.reacs
    for reaction, flux in zip(reactions, fluxResults.flux):
        if not reaction.name.endswith("_transp"):
            changeDic[reaction.name] = float('%.3g' % flux)

    vList = changeDic.values()

    for i in xrange(len(vList)):
        vList[i] = math.sqrt(math.pow(vList[i], 2))
    vList = sorted(vList)

    oldMin = vList[0]
    oldMax = vList[-1]
    oldRange = oldMax - oldMin
    newMin = 1
    newMax = 10
    newRange = newMax - newMin

    enzymes = org.enzymes
    nodeDic = OrderedDict()
    nodecounter = 0
    graph = nx.DiGraph()
    for enzyme in enzymes:
        if not enzyme.name.endswith("_transp"):
            nodecounter += 1
            nodeDic[enzyme.name] = nodecounter

            for substrate in enzyme.substrates:
                nodecounter += 1
                if substrate not in nodeDic:
                    nodeDic[substrate] = nodecounter
                    graph.add_node(nodeDic[substrate], label=substrate, shape="oval", color=str((nodecounter % 8)+1))

            for product in enzyme.products:
                nodecounter += 1
                if product not in nodeDic:
                    nodeDic[product] = nodecounter
                    graph.add_node(nodeDic[product], label=product, shape="oval", color=str((nodecounter % 8)+1))

            for substrate in enzyme.substrates:
                for product in enzyme.products:
                    flux = changeDic[enzyme.name]

                    color = "black"
                    if flux < 0:
                        color = "red"
                    elif flux > 0:
                        color = "green"

                    value = flux
                    if value != 0:
                        value = (math.sqrt(math.pow(value, 2)) - oldMin) / oldRange * newRange + newMin

                    attr = {"color": color, "penwidth": value, "label": "{} ({})".format(enzyme.name, flux), "object": enzyme}

                    graph.add_edge(nodeDic[substrate], nodeDic[product], attr)
                    if enzyme.reversible:
                        graph.add_edge(nodeDic[product], nodeDic[substrate], label=enzyme.name, object=enzyme)

    return graph, nodeDic
