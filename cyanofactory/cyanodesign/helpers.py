"""
Copyright (c) 2014 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""
from collections import OrderedDict
import math
from networkx.readwrite import json_graph
from PyNetMet2.enzyme import Enzyme
from PyNetMet2.metabolism import Metabolism
from networkx import DiGraph

def apply_commandlist(model, commandlist):
    """
    :type model: PyNetMet2.metabolism.Metabolism
    """
    for command in commandlist:
        try:
            op = command["op"]
            typ = command["type"]
            name = command["name"]
            obj = command["object"]
        except KeyError:
            raise ValueError("Bad command " + str(command))

        if typ == "reaction":
            enzyme = model.get_reaction(name)

            substrates = []
            products = []

            if op in ["add", "edit"]:
                if "name" not in obj:
                    raise ValueError("Bad commmand " + str(command))

                try:
                    for m, s in zip(obj["metabolites"], obj["stoichiometry"]):
                        (substrates if s < 0 else products).append([str(abs(s)), m])
                except KeyError:
                    pass

            if op == "add":
                if enzyme is not None:
                    raise ValueError("Reaction already in model: " + name)
                if name != obj["name"]:
                    raise ValueError("Reaction name mismatch: " + name)
                if "reversible" not in obj:
                    raise ValueError("Bad command " + str(command))

                reac = "{} : {} {} {}".format(
                    obj["name"],
                    " + ".join([" ".join(x) for x in substrates]),
                    "<->" if obj["reversible"] else "->",
                    " + ".join([" ".join(x) for x in products])
                )

                model.add_reaction(reac)
                reaction = model.get_reaction(name)

                try:
                    reaction.constraint = obj["constraints"]
                except KeyError:
                    pass

                try:
                    reaction.pathway = obj["pathway"]
                except KeyError:
                    pass

                try:
                    reaction.disabled = not obj["enabled"]
                except KeyError:
                    pass

                model.calcs()
            elif op == "edit":
                if enzyme is None:
                    raise ValueError("Reaction not in model: " + name)

                if "metabolites" in obj and "substrates" in obj:
                    enzyme.substrates = map(lambda x: x[1], substrates)
                    enzyme.products = map(lambda x: x[1], products)
                    enzyme.stoic[0] = map(lambda x: x[0], substrates)
                    enzyme.stoic[1] = map(lambda x: x[1], products)

                try:
                    enzyme.reversible = obj["reversible"]
                except KeyError:
                    pass

                try:
                    enzyme.constraint = obj["constraints"]
                except KeyError:
                    pass

                try:
                    enzyme.pathway = obj["pathway"]
                except KeyError:
                    pass

                try:
                    enzyme.disabled = not obj["enabled"]
                except KeyError:
                    pass

                if obj["name"] != enzyme.name and model.has_reaction(obj["name"]):
                    raise ValueError("Reaction already in model: " + obj["name"])

                enzyme.name = obj["name"]

                model.calcs()
            elif op == "delete":
                if enzyme is None:
                    raise ValueError("Reaction not in model: " + name)

                model.remove_reaction(enzyme.name)
            else:
                raise ValueError("Invalid operation " + op)

        elif typ == "metabolite":
            if "name" not in obj:
                raise ValueError("Bad commmand " + str(command))

            if op in ["add", "edit", "delete"]:
                model.rename_metabolite(name, obj["name"])

                if "external" in obj:
                    if obj["external"] and op != "delete":
                        model.make_metabolite_external(obj["name"])
                    else:
                        # These are deleted when non references them
                        model.make_metabolite_internal(obj["name"])
            else:
                raise ValueError("Invalid operation " + op)

        elif typ == "pathway":
            if "name" not in obj:
                raise ValueError("Bad commmand " + str(command))

            if op == "edit":
                for reac in model.enzymes:
                    if reac.pathway == name:
                        reac.pathway = obj["name"]
                model.calcs()
            else:
                raise ValueError("Invalid operation " + op)

        else:
            raise ValueError("Invalid command " + typ)

    model.calcs()
    return model

def model_from_string(model_str):
    import os
    import tempfile

    with tempfile.NamedTemporaryFile(delete=False) as fid:
        path = fid.name

        fid.write(model_str.encode("utf-8"))

    try:
        org = Metabolism(path)
    except:
        return None
    finally:
        os.remove(path)

    return org

def flatten(li):
    return [item for sublist in li for item in sublist]

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
        metabolites += org.get_reaction(reac).metabolites

    met_ids = map(lambda x: nodeDic[x], metabolites)
    g = json_graph.node_link_graph(jsonGraph)

    g.remove_edges_from(filter(lambda x: g.get_edge_data(*x)["object"].name not in reacIDs, g.edges(met_ids)))

    # Get products/substrates directly connected to filter
    #reacIDs += flatten(g.in_edges(reacIDs)) + flatten(g.out_edges(reacIDs))

    h = g.subgraph(met_ids)
    return h

def calc_reactions(org, fluxResults):
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

    enzymes = fluxResults.reacs
    nodeDic = OrderedDict()
    nodecounter = 0
    graph = DiGraph()

    objective = None
    if len(org.obj) > 0:
        objective = org.obj[0][0]

    for enzyme in enzymes:
        if not enzyme.name.endswith("_transp"):
            nodecounter += 1
            nodeDic[enzyme.name] = nodecounter
            if enzyme.name == objective:
                graph.add_node(nodeDic[enzyme.name])

            for substrate in enzyme.substrates:
                nodecounter += 1
                if substrate not in nodeDic:
                    nodeDic[substrate] = nodecounter
                    attr = {"label": substrate, "shape": "oval", "color": str((nodecounter % 8)+1)}
                    graph.add_node(nodeDic[substrate], **attr)

                if enzyme.name == objective:
                    graph.add_nodes_from([(nodeDic[substrate], {"shape":"box"})])

            for product in enzyme.products:
                nodecounter += 1
                if product not in nodeDic:
                    nodeDic[product] = nodecounter
                    attr = {"label": product, "shape": "oval", "color": str((nodecounter % 8)+1)}
                    graph.add_node(nodeDic[product], **attr)

                if enzyme.name == objective:
                    graph.add_nodes_from([(nodeDic[product], {"shape":"box"})])

            # Check if enzyme is objective
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
                    else:
                        value = 1

                    attr = {"color": color, "penwidth": value, "label": "{} ({})".format(enzyme.name, flux), "object": enzyme}

                    if enzyme.name == objective:
                        attr["style"] = "dashed"

                    graph.add_edge(nodeDic[substrate], nodeDic[product], attr)

                    if enzyme.name == objective:
                        graph.add_edge(nodeDic[substrate], nodeDic[enzyme.name])
                        graph.add_edge(nodeDic[enzyme.name], nodeDic[product])

                    if enzyme.reversible:
                        graph.add_edge(nodeDic[product], nodeDic[substrate], label=enzyme.name, object=enzyme)

    return graph, nodeDic
