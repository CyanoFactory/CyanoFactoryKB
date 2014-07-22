import json
from django.contrib.auth.decorators import login_required
from django.http.response import HttpResponse, HttpResponseBadRequest
from cyano.decorators import ajax_required
from PyNetMet.metabolism import *
from PyNetMet.fba import *
from django.conf import settings
from cyano.helpers import render_queryset_to_response
import networkx as nx
from networkx.readwrite import json_graph
import re
import math
import pygraphviz
import json


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

    graphInfos = drawDesign(org)
    graph = nx.to_agraph(json_graph.node_link_graph(graphInfos[0]))
    graph.graph_attr.update(splines=True, overlap=False, rankdir="LR")
    graph.node_attr.update(style="filled", colorscheme="pastel19")
    outgraph = str(graph.to_string())
    return HttpResponse(
        json.dumps({"external": org.external, "enzymes": ret, "objective": org.objective, "graph": outgraph}),
        content_type="application/json")


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

  # print "\n".join(data)
  # print constr
  # print ext
  # print objective

    try:
        organism = Metabolism("model_name",
                              reactions=data,
                              constraints=constr,
                              external=ext,
                              objective=objective,
                              fromfile=False)

        fba = FBA(organism)

        graphInfo = drawDesign(organism)
        xgraph = graphInfo[0]
        nodeIDs = graphInfo[1]
        g = calcReactions(xgraph, nodeIDs, str(fba))
        graph = nx.to_agraph(g)
        graph.graph_attr.update(splines=True, overlap=False, rankdir="LR")
        graph.node_attr.update(style="filled", colorscheme="pastel19")
        outgraph = str(graph.to_string())

    except ValueError:
        return HttpResponseBadRequest("Invalid data type")

    return HttpResponse(outgraph)


def export(request):
    pass


def drawDesign(org):
    """
    Generating Graph for all Reaction from defined organism.
    return networkx Graph in JSON-Format and Node-ID - Node-Label Dictionary[Node-Label: Node-ID]
    @param org: object organism from PyNetMet
    @return: Dictionary["graph": Graph in JSON-Format, "nodeDic": Dictionary of Node-IDs]
    """
    enzymes = org.enzymes
    nodeDic = {}
    nodecounter = 0
    graph = nx.DiGraph()
    for enzyme in enzymes:
        if not enzyme.name.endswith("_transp"):
            nodecounter += 1
            nodeDic[enzyme.name] = nodecounter
            graph.add_node(nodeDic[enzyme.name], label=enzyme.name, shape="box")
            for substrate in enzyme.substrates:
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


def getSelectedReaction(reacIDs, jsonGraph):
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
    @return Graph with added Attributes
    """
    print nodeDic
    changeDic = readResults(fluxResults)
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

def readResults(resultText):
    """
    Generating a Dictionary from the PyNetMet FLUX-Analysis with a given String Input.
    @param resultText: Output String from FLUX-Analysis
    @return: Dictionary[Substrate: FLUX-Result]
    """
    searchP = r"(\S+)\s+---->\s+(-?\d+\.\d+)"
    #print resultText
    resultList = re.findall(searchP, resultText)
    changeDic = {}
    for result in resultList:
        if not result[0].endswith("_transp"):
            changeDic[result[0]] = float(result[1])

    return changeDic

