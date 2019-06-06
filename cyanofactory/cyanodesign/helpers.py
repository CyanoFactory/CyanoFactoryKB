"""
Copyright (c) 2014 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""
from collections import OrderedDict
import math
from networkx.readwrite import json_graph
from networkx import DiGraph
from metabolic_model import metabolic_model

def get_selected_reaction(jsonGraph, nodeDic, reacIDs, org):
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
        r = org.reaction.get(id=reac)
        metabolites += r.substrates
        metabolites += r.products

    met_ids = list(map(lambda x: nodeDic[x.id], metabolites))
    g = json_graph.node_link_graph(jsonGraph)

    def l(x):
        return g.get_edge_data(*x)["object"]["id"] not in reacIDs

    g.remove_edges_from(list(filter(l, g.edges(met_ids))))

    # Get products/substrates directly connected to filter
    #reacIDs += flatten(g.in_edges(reacIDs)) + flatten(g.out_edges(reacIDs))

    h = g.subgraph(met_ids)
    return h

def calc_reactions(org: metabolic_model.MetabolicModel, fluxResults: metabolic_model.SimulationResult, objective):
    """
    Adding Results from FLUX-Analysis
    @param org: organism
    @param fluxResults: Results from PyNetMet FLUX-Analysis Calculation as String
    @return Graph with added Attributes
    """
    changeDic = {}
    for result in fluxResults.results:
        if not result.reaction.endswith("_transp"):
            changeDic[result.reaction] = float('%.3g' % result.flux)

    vList = list(changeDic.values())

    for i in range(len(vList)):
        vList[i] = math.sqrt(math.pow(vList[i], 2))
    vList = sorted(vList)

    oldMin = vList[0]
    oldMax = vList[-1]
    oldRange = oldMax - oldMin
    newMin = 1
    newMax = 10
    newRange = newMax - newMin

    enzymes = org.reactions
    nodeDic = OrderedDict()
    nodecounter = 0
    graph = DiGraph()

    for enzyme in enzymes:
        if enzyme.enabled and not enzyme.id.endswith("_transp"):
            nodecounter += 1
            nodeDic[enzyme.id] = nodecounter
            if enzyme == objective:
                graph.add_node(nodeDic[enzyme.id])

            for substrate in enzyme.substrates:
                nodecounter += 1
                if substrate.id not in nodeDic:
                    nodeDic[substrate.id] = nodecounter
                    attr = {"label": org.metabolite.get(id=substrate.id).name, "shape": "oval", "color": str((nodecounter % 8)+1)}
                    graph.add_node(nodeDic[substrate.id], **attr)

                if enzyme.id == objective:
                    graph.add_nodes_from([(nodeDic[substrate.id], {"shape":"box"})])

            for product in enzyme.products:
                nodecounter += 1
                if product.id not in nodeDic:
                    nodeDic[product.id] = nodecounter
                    attr = {"label": org.metabolite.get(id=product.id).name, "shape": "oval", "color": str((nodecounter % 8)+1)}
                    graph.add_node(nodeDic[product.id], **attr)

                if enzyme.id == objective:
                    graph.add_nodes_from([(nodeDic[product.id], {"shape":"box"})])

            # Check if enzyme is objective
            for substrate in enzyme.substrates:
                for product in enzyme.products:
                    flux = changeDic[enzyme.id]

                    color = "black"
                    if flux < 0:
                        color = "red"
                    elif flux > 0:
                        color = "green"

                    value = flux
                    if value != 0:
                        value = (math.sqrt(math.pow(value, 2)) - oldMin) / oldRange * newRange + newMin
                    else:
                        value = 1

                    attr = {"color": color, "penwidth": value, "label": u"{} ({})".format(enzyme.name, flux), "object": {"name":enzyme.name,"id":enzyme.id}}

                    if enzyme == objective:
                        attr["style"] = "dashed"

                    graph.add_edge(nodeDic[substrate.id], nodeDic[product.id], attr)

                    if enzyme == objective:
                        graph.add_edge(nodeDic[substrate.id], nodeDic[enzyme.id])
                        graph.add_edge(nodeDic[enzyme.id], nodeDic[product.id])

                    if enzyme.reversible:
                        graph.add_edge(nodeDic[product.id], nodeDic[substrate.id], label=enzyme.id, object={"name":enzyme.name,"id":enzyme.id})

    return graph, nodeDic
