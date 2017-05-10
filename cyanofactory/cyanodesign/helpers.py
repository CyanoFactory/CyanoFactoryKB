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
from cyanodesign.json_model import JsonModel
from PyNetMet2.util import create_sid

def apply_commandlist(model, commandlist):
    """
    :type model: PyNetMet2.Metabolism
    """
    for command in commandlist:
        try:
            op = command["op"]
            typ = command["type"]
            name = command["id"]
            cid = create_sid(name)
            obj = command["object"]
        except KeyError:
            raise ValueError("Bad command " + str(command))

        if typ == "reaction":
            enzyme = model.get_reaction(name)

            if op in ["add", "edit"]:
                if "name" not in obj:
                    raise ValueError("Bad commmand " + str(command))

            if op == "add":
                if enzyme is not None:
                    raise ValueError("Reaction already in model: " + name)
                if name != obj["id"]:
                    raise ValueError("Reaction name mismatch: " + name)
                if "reversible" not in obj:
                    raise ValueError("Bad command " + str(command))

                model.add_reaction(JsonModel.Reaction(**obj))
            elif op == "edit":
                if enzyme is None:
                    raise ValueError("Reaction not in model: " + name)

                if obj["id"] != enzyme.id and model.has_reaction(obj["id"]):
                    raise ValueError("Reaction already in model: " + obj["name"])

                enzyme.id = obj["id"]
                enzyme.name = obj["name"]

                if "substrates" in obj:
                    enzyme.substrates = map(lambda x: JsonModel.Compound(**x), obj["substrates"])
                if "products" in obj:
                    enzyme.products = map(lambda x: JsonModel.Compound(**x), obj["products"])

                try:
                    enzyme.reversible = obj["reversible"]
                except KeyError:
                    pass

                try:
                    enzyme.constraints = obj["constraints"]
                except KeyError:
                    pass

                try:
                    enzyme.pathway = obj["pathway"]
                except KeyError:
                    pass

                try:
                    enzyme.disabled = obj["disabled"]
                except KeyError:
                    pass
            elif op == "delete":
                if enzyme is None:
                    raise ValueError("Reaction not in model: " + name)

                model.remove_reaction(enzyme.name)
            else:
                raise ValueError("Invalid operation " + op)

        elif typ == "metabolite":
            if op in ["add", "edit"] and "id" not in obj:
                raise ValueError("Bad command " + str(command))

            if op == "add":
                if name != obj["id"]:
                    raise ValueError("Metabolite name mismatch: " + name)

                model.add_metabolite(obj["id"], obj.get("external", False))
            elif op == "edit":
                if "external" in obj:
                    if obj["external"]:
                        model.make_metabolite_external(name)
                    else:
                        model.make_metabolite_internal(name)
                if name != obj["name"]:
                    model.rename_metabolite(name, obj["name"])
            elif op == "delete":
                model.remove_metabolite(name)
            else:
                raise ValueError("Invalid operation " + op)

        elif typ == "pathway":
            if "name" not in obj:
                raise ValueError("Bad commmand " + str(command))

            if op == "edit":
                for reac in model.enzymes:
                    if reac.pathway == name:
                        reac.pathway = obj["id"]

            else:
                raise ValueError("Invalid operation " + op)

        else:
            raise ValueError("Invalid command " + typ)

    return model

def compress_command_list(commandlist):
    """
    :type model: json_model.JsonModel
    """
    last_command = None

    new_commandlist = []
    pending_commands = []

    # No error checking use after apply_commandlist

    for command in commandlist:
        op = command["op"]
        typ = command["type"]
        name = command["id"]

        if len(pending_commands) == 0:
            pending_commands.append(command)
            last_command = command
            continue

        if typ == "metabolite" or typ == "reaction" or typ == "pathway":
            if last_command["type"] == typ:
                if last_command["object"].get("id", last_command["id"]) != name or op == "add":
                    # Different to previous ones
                    new_commandlist += pending_commands
                    pending_commands = []
                    pending_commands.append(command)
                elif op == "edit":
                    if last_command["object"].get("id", last_command["id"]) == name:
                        # Merge with previous
                        pending_commands[-1]["object"].update(command["object"])
                        # If previous was an add update the name
                        if pending_commands[-1]["op"] == "add":
                            pending_commands[-1]["id"] = command["object"]["id"]
                    else:
                        new_commandlist += pending_commands
                        pending_commands = []
                elif op == "delete":
                    # Remove all previous instances and take name of first
                    command["id"] = pending_commands[0]["id"]
                    # If first command was an "add" the item was created and just deleted
                    # can be removed completely
                    if pending_commands[0]["op"] != "add":
                        new_commandlist.append(command)

                    pending_commands = []
            else:
                new_commandlist += pending_commands
                pending_commands = []
                pending_commands.append(command)

        last_command = command

    new_commandlist += pending_commands
    return new_commandlist

def model_from_string(model_str):
    return JsonModel.from_json(model_str)

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

    met_ids = list(map(lambda x: nodeDic[x], metabolites))
    g = json_graph.node_link_graph(jsonGraph)

    g.remove_edges_from(list(filter(lambda x: g.get_edge_data(*x)["object"]["name"] not in reacIDs, g.edges(met_ids))))

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

                    attr = {"color": color, "penwidth": value, "label": u"{} ({})".format(enzyme.name, flux), "object": {"name":enzyme.name,"id":enzyme.id}}

                    if enzyme.name == objective:
                        attr["style"] = "dashed"

                    graph.add_edge(nodeDic[substrate], nodeDic[product], attr)

                    if enzyme.name == objective:
                        graph.add_edge(nodeDic[substrate], nodeDic[enzyme.name])
                        graph.add_edge(nodeDic[enzyme.name], nodeDic[product])

                    if enzyme.reversible:
                        graph.add_edge(nodeDic[product], nodeDic[substrate], label=enzyme.name, object={"name":enzyme.name,"id":enzyme.id})

    return graph, nodeDic
