"""
Copyright (c) 2014 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

import StringIO
from collections import OrderedDict
import json
import os
import tempfile
from django.contrib.auth.decorators import login_required
from django.core.exceptions import ObjectDoesNotExist
from django.http.response import HttpResponse, HttpResponseBadRequest
from django.shortcuts import redirect
from django.views.decorators.csrf import ensure_csrf_cookie
from networkx.algorithms.shortest_paths.generic import has_path, shortest_path, all_shortest_paths
from networkx.exception import NetworkXNoPath
from cyano.decorators import ajax_required
from PyNetMet.metabolism import *
from PyNetMet.fba import *
from cyano.helpers import render_queryset_to_response, render_queryset_to_response_error
from cyanodesign.forms import UploadModelForm
from cyanodesign.helpers import model_from_string, apply_commandlist
from .models import DesignModel
import networkx as nx
from networkx.readwrite import json_graph
import re
import math
import pygraphviz
import json

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

    for enzyme in org.reactions:
        # Filter out auto-created external transport funcs
        if not enzyme.name.endswith("_ext_transp"):
            ret.append({
                "name": enzyme.name,
                "stoichiometric": [enzyme.reactants_stoic] + [enzyme.products_stoic],
                "substrates": map(lambda m: m.name, enzyme.reactants),
                "products": map(lambda m: m.name, enzyme.products),
                "reversible": enzyme.reversible,
                "constraints": enzyme.constraint,
                "disabled": enzyme.disabled
            })

    #graphInfos = drawDesign(org)
    #graph = nx.to_agraph(json_graph.node_link_graph(graphInfos[0]))
    #graph.graph_attr.update(splines=True, overlap=False, rankdir="LR")
    #graph.node_attr.update(style="filled", colorscheme="pastel19")
    #outgraph = str(graph.to_string())
    outgraph = ""

    return HttpResponse(json.dumps({
        "external": map(lambda m: m.name, filter(lambda m: m.external, org.get_metabolites())),
        "enzymes": ret,
        "objective": org.objective.name if org.objective else None,
        "graph": outgraph}), content_type="application/json")

@login_required
@ajax_required
def simulate(request, pk):
    if not all(x in request.POST for x in ["commands", "disabled", "objective", "display", "auto_flux"]):
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
        obj_reac = org.get_reaction(objective)
        if obj_reac is None:
            raise ValueError("Objective not in model: " + objective)
        org.objective = obj_reac
    except ValueError as e:
        return HttpResponseBadRequest("Model error: " + e.message)

    for reaction in org.reactions:
        reaction.disabled = False

    for item in disabled:
        try:
            reac = org.get_reaction(item)
        except ValueError:
            return HttpResponseBadRequest("Bad disabled list: " + item)

        reac.disabled = True

    if obj_reac.disabled:
        return HttpResponseBadRequest("Objective disabled: " + objective)

    for item in display:
        if not org.has_reaction(item):
            return HttpResponseBadRequest("Bad display list: " + item)

    try:
        org.fba()
    except ValueError as e:
        return HttpResponseBadRequest("FBA error: " + e.message)

    graphInfo = drawDesign(org)
    xgraph = graphInfo[0]
    nodeIDs = graphInfo[1]
    full_g = calcReactions(xgraph, nodeIDs, org)

    import itertools

    display = json.loads(request.POST["display"])

    # Auto filter by FBA
    if auto_flux:
        flux = filter(lambda x: not x.disabled, org.reactions[:])
        flux = sorted(flux, key=lambda x: -x.flux)
        display = map(lambda x: x.name, flux[:10])

    g = getSelectedReaction(json_graph.node_link_data(full_g), nodeIDs, display)

    # Add Pseudopaths
    for reac_pairs in itertools.combinations(display, 2):
        if not has_path(g, nodeIDs[reac_pairs[0]], nodeIDs[reac_pairs[1]]):
            # Test for theoretical path
            try:
                all_paths = all_shortest_paths(full_g, nodeIDs[reac_pairs[0]], nodeIDs[reac_pairs[1]])

                for path in all_paths:
                    # Take second to last
                    #g.add_edge(path[1], path[-2])
                    #print "Pseudopath from", path[1], "to", path[-2]
                    pass
            except NetworkXNoPath:
                pass

        if not has_path(g, nodeIDs[reac_pairs[1]], nodeIDs[reac_pairs[0]]):
            # Test for theoretical path
            try:
                all_paths = all_shortest_paths(full_g, nodeIDs[reac_pairs[1]], nodeIDs[reac_pairs[0]])

                for path in all_paths:
                    # Take second to last
                    #g.add_edge(path[1], path[-2])
                    #print "Pseudopath from", path[1], "to", path[-2]
                    pass
            except NetworkXNoPath:
                pass

    graph = nx.to_agraph(g)
    graph.graph_attr.update(splines=True, overlap=False, rankdir="LR")
    graph.node_attr.update(style="filled", colorscheme="pastel19")
    outgraph = str(graph.to_string())

    return HttpResponse(json.dumps(
        {"graph": outgraph,
        "fluxes": map(lambda x: [x.name, x.flux], org.reactions)}),
        content_type="application/json"
    )


@login_required
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


@login_required
@ajax_required
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
        obj_reac = org.get_reaction(objective)
        if obj_reac is None:
            raise ValueError("Objective not in model: " + objective)
        org.objective = obj_reac
    except ValueError as e:
        return HttpResponseBadRequest("Model error: " + e.message)

    for reaction in org.reactions:
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


def drawDesign(org):
    """
    Generating Graph for all Reaction from defined organism.
    return networkx Graph in JSON-Format and Node-ID - Node-Label Dictionary[Node-Label: Node-ID]
    @param org: object organism from PyNetMet
    @return: Dictionary["graph": Graph in JSON-Format, "nodeDic": Dictionary of Node-IDs]
    """
    enzymes = org.reactions
    nodeDic = OrderedDict()
    nodecounter = 0
    graph = nx.DiGraph()
    for enzyme in enzymes:
        if not enzyme.name.endswith("_transp"):
            nodecounter += 1
            nodeDic[enzyme.name] = nodecounter
            graph.add_node(nodeDic[enzyme.name], label=enzyme.name, shape="box")
            for substrate in enzyme.reactants:
                substrate = substrate.name
                nodecounter += 1
                if not substrate in nodeDic:
                    nodeDic[substrate] = nodecounter
                    graph.add_node(nodeDic[substrate], label=substrate, shape="oval", color=str((nodecounter % 8)+1))
                graph.add_edge(nodeDic[substrate], nodeDic[enzyme.name])
                if enzyme.reversible:
                    graph.add_edge(nodeDic[enzyme.name], nodeDic[substrate])
            for product in enzyme.products:
                product = product.name
                nodecounter += 1
                if not product in nodeDic:
                    nodeDic[product] = nodecounter
                    graph.add_node(nodeDic[product], label=product, shape="oval", color=str((nodecounter % 8)+1))
                graph.add_edge(nodeDic[enzyme.name], nodeDic[product])
                if enzyme.reversible:
                    graph.add_edge(nodeDic[enzyme.name], nodeDic[product])

    graphJson = json_graph.node_link_data(graph)
    data = [graphJson, nodeDic]

    return data

def flatten(li):
    return [item for sublist in li for item in sublist]

def getSelectedReaction(jsonGraph, nodeDic, reacIDs):
    """
    Filtering selected Reactions and show Results from PyNetMet calculation.
    It returns a subgraph of the Graph from the jsonGraph. The output is a DOT-Language String.
    @param jsonGraph: Graph in JSON-Format
    @param nodeDic: dict mapping names to ids
    @param reacIDs: Name of reactions that contained in the nodeDic
    @return Subgraph
    """
    # Translate reac names to IDs
    reacIDs = map(lambda x: nodeDic[x], reacIDs)
    g = json_graph.node_link_graph(jsonGraph)

    # Get products/substrates directly connected to filter
    reacIDs += flatten(g.in_edges(reacIDs)) + flatten(g.out_edges(reacIDs))

    h = g.subgraph(reacIDs)
    return h

def calcReactions(jsonGraph, nodeDic, fluxResults):
    """
    Adding Results from FLUX-Analysis
    @param jsonGraph: Graph in JSON-Format
    @param nodeDic: Dictionary with keys = Node-Name and value = node_id
    @param fluxResults: Results from PyNetMet FLUX-Analysis Calculation as String
    @type fluxResults: bioparser.optgene.OptGeneParser
    @return Graph with added Attributes
    """
    changeDic = {}
    reactions = filter(lambda x: not x.disabled, fluxResults.reactions)
    for reaction in reactions:
        if not reaction.name.endswith("_transp"):
            changeDic[reaction.name] = float('%.3g' % reaction.flux)

    g = json_graph.node_link_graph(jsonGraph)
    vList = changeDic.values()

    for i in xrange(len(vList)):
        vList[i] = math.sqrt(math.pow(vList[i], 2))
    vList = sorted(vList)

    if len(vList) == 0:
        return g

    oldMin = vList[0]
    oldMax = vList[-1]
    oldRange = oldMax - oldMin
    newMin = 1
    newMax = 10
    newRange = newMax - newMin

    for key in changeDic:
        node = nodeDic[key]
        value = changeDic[key]

        color = "black"
        if value < 0:
            color = "red"
        elif value > 0:
            color = "green"

        newValue = 1
        if value != 0:
            newValue = (math.sqrt(math.pow(value, 2)) - oldMin) / oldRange * newRange + newMin

        thikness = newValue

        edgeList = g.out_edges(node)
        edgeList.extend(g.in_edges(node))
        for aEdge in edgeList:
            g.edge[aEdge[0]][aEdge[1]]["color"] = color
            g.edge[aEdge[0]][aEdge[1]]["penwidth"] = thikness
            g.edge[aEdge[0]][aEdge[1]]["label"] = value

    return g
