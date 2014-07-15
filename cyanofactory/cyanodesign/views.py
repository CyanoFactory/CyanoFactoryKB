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

    graph = drawDesign(org)
    return HttpResponse(
        json.dumps({"external": org.external, "enzymes": ret, "objective": org.objective, "graph": graph}),
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

    print "\n".join(data)
    print constr
    print ext
    print objective

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


def drawDesign(org):
    enzymes = org.enzymes
    svg_graph = "digraph g{k=10;\n splines=true;\n overlap=false;\n edge [dir = both];"
    nodeDic = {}
    edgeDic = {}
    nodecounter = 0
    edgecounter = 0
    graph = nx.Graph()
    svg_graph += "node [colorscheme = pastel19, style = filled];\n"
    for enzyme in enzymes:
        # TODO better Selection of not usable Nodes
        if not "_transp" in enzyme.name:
            nodecounter += 1
            nodeDic[enzyme.name] = nodecounter
            tail = "none"
            graph.add_node(nodeDic[enzyme.name], label=enzyme.name, enzyme=True, reversible=enzyme.reversible)
            svg_graph += str(nodeDic[enzyme.name]) + '[label = "' + enzyme.name + '", shape = box, color=' + str(
                (nodecounter % 8) + 1) + '];\n'
            if enzyme.reversible:
                tail = "normal"
            for substrate in enzyme.substrates:
                nodecounter += 1
                if not substrate in nodeDic:
                    nodeDic[substrate] = nodecounter
                    graph.add_node(nodeDic[substrate], label=substrate, enzyme=False)
                    svg_graph += str(nodeDic[substrate]) + '[label = "' + substrate + '"];\n'
                graph.add_edge(nodeDic[substrate], nodeDic[enzyme.name])
                if enzyme.reversible:
                    graph.add_edge(nodeDic[enzyme.name], nodeDic[substrate])
                svg_graph = svg_graph + " " + str(nodeDic[substrate]) + " -> " + str(
                    nodeDic[enzyme.name]) + "[arrowtail=" + tail + ", arrowhead=normal];\n"
            for product in enzyme.products:
                nodecounter += 1
                if not product in nodeDic:
                    nodeDic[product] = nodecounter
                    graph.add_node(nodeDic[product], label=product, enzyme=False)
                    svg_graph += str(nodeDic[product]) + '[label = "' + product + '"];\n'
                graph.add_edge(nodeDic[product], nodeDic[enzyme.name])
                if enzyme.reversible:
                    graph.add_edge(nodeDic[enzyme.name], nodeDic[product])
                svg_graph = svg_graph + " " + str(nodeDic[enzyme.name]) + " -> " + str(
                    nodeDic[product]) + "[arrowtail=" + tail + ", arrowhead=normal];\n"
    svg_graph += "}"
    #print svg_graph
    print json_graph.node_link_data(graph)
    return svg_graph

