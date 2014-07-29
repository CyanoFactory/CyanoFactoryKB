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
                "constraints": enzyme.constraint
            })

    graphInfos = drawDesign(org)
    graph = nx.to_agraph(json_graph.node_link_graph(graphInfos[0]))
    graph.graph_attr.update(splines=True, overlap=False, rankdir="LR")
    graph.node_attr.update(style="filled", colorscheme="pastel19")
    outgraph = str(graph.to_string())

    return HttpResponse(json.dumps({
        "external": map(lambda m: m.name, filter(lambda m: m.external, org.get_metabolites())),
        "enzymes": ret,
        "objective": org.objective.name if org.objective else None,
        "graph": outgraph}), content_type="application/json")

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

    if not all(isinstance(x, list) for x in [commands, disabled]):
        return HttpResponseBadRequest("Invalid data type")

    if not isinstance(objective, basestring):
        return HttpResponseBadRequest("Invalid data type")

    try:
        org = apply_commandlist(org, commands)
        obj_reac = org.get_reaction(objective)
        if obj_reac is None:
            raise ValueError("Objective not in model: " + objective)
    except ValueError as e:
        return HttpResponseBadRequest("Model error: " + e.message)

    try:
        org.fba()
    except ValueError as e:
        return HttpResponseBadRequest("FBA error: " + e.message)

    graphInfo = drawDesign(org)
    xgraph = graphInfo[0]
    nodeIDs = graphInfo[1]
    g = calcReactions(xgraph, nodeIDs, org)
    graph = nx.to_agraph(g)
    graph.graph_attr.update(splines=True, overlap=False, rankdir="LR")
    graph.node_attr.update(style="filled", colorscheme="pastel19")
    outgraph = str(graph.to_string())

    return HttpResponse(outgraph)



    #return HttpResponse(json.dumps({"flux": map(lambda x: [x.name, x.flux], org.reactions)}), content_type="application/json")

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


def drawDesign(org):
    """
    Generating Graph for all Reaction from defined organism.
    return networkx Graph in JSON-Format and Node-ID - Node-Label Dictionary[Node-Label: Node-ID]
    @param org: object organism from PyNetMet
    @return: Dictionary["graph": Graph in JSON-Format, "nodeDic": Dictionary of Node-IDs]
    """
    enzymes = org.reactions
    nodeDic = {}
    nodecounter = 0
    graph = nx.DiGraph()
    for enzyme in enzymes:
        if not enzyme.name.endswith("_transp"):
            nodecounter += 1
            nodeDic[enzyme.name] = nodecounter
            graph.add_node(nodeDic[enzyme.name], label=enzyme.name, shape="box")
            for substrate in enzyme.reactants:
                nodecounter += 1
                if not substrate in nodeDic:
                    nodeDic[substrate] = nodecounter
                    graph.add_node(nodeDic[substrate], label=substrate, shape="oval", color=str((nodecounter % 8)+1))
                graph.add_edge(nodeDic[substrate], nodeDic[enzyme.name])
                if enzyme.reversible:
                    graph.add_edge(nodeDic[enzyme.name], nodeDic[substrate])
            for product in enzyme.products:
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


def getSelectedReaction(reacIDs, nodeDic, jsonGraph):
    """
    Filtering selected Reactions and show Results from PyNetMet calculation.
    It returns a subgraph of the Graph from the jsonGraph. The output is a DOT-Language String.
    @param reacIDs: Array of Reaction IDs which are contained in the Graph of the jsonFile
    @param jsonGraph: Graph in JSON-Format
    @return Subgraph in DOT-Language
    """

    g = json_graph.node_link_graph(jsonGraph)
    prelistofNodes = []
    for reacID in reacIDs:
        prelistofNodes = g[reacID]

    listofNodes = []
    for node in prelistofNodes:
        listofNodes.extend(g[node])
    listofNodes.extend(reacIDs)

    h = g.subgraph(listofNodes)
    hGraphviz = str(nx.to_agraph(h).to_string())
    return hGraphviz


def calcReactions(jsonGraph, nodeDic, fluxResults):
    """
    Adding Results from FLUX-Analysis
    @param jsonGraph: Graph in JSON-Format
    @param nodeDic: Dictionary with keys = Node-Name and value = node_id
    @param fluxResults: Results from PyNetMet FLUX-Analysis Calculation as String
    @type fluxResults: bioparser.optgene.OptGeneParser
    @return Graph with added Attributes
    """
    print nodeDic
    changeDic = {}
    reactions = fluxResults.reactions
    for reaction in reactions:
        if not reaction.name.endswith("_transp"):
            changeDic[reaction.name] = reaction.flux

    g = json_graph.node_link_graph(jsonGraph)
    vList = changeDic.values()
    while 0 in vList:
        vList.remove(0)

    for i in xrange(len(vList)):
        vList[i] = math.sqrt(math.pow(vList[i], 2))
    vList = sorted(vList)

    oldMin = vList[0]
    oldMax = vList[-1]
    oldRange = oldMax - oldMin
    newMin = 1
    newMax = 20
    newRange = newMax - newMin

    for key in changeDic:
        node = nodeDic[key]
        value = changeDic[key]

        print value
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

