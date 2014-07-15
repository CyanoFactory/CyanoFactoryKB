from django.core.exceptions import ObjectDoesNotExist
from django.http import HttpResponse
from django.shortcuts import render_to_response
from cyano.decorators import ajax_required
from stringdb.models import Proteins, NodeNodeLinks, ProteinsOrthgroups, Chemicals, ProteinChemicalLinks, \
    ChemicalChemicalLinks, ProteinChemicalLinksDetailed, ChemicalChemicalLinksDetailed
import networkx as nx
import json
from networkx.readwrite import json_graph
import os
from cyano.models import ProteinComparison
from cyano.helpers import render_queryset_to_response
import re

''' Getting Interaction Information for selected ProteinID '''


def checkInteraction(request, protID, limit=10, chemlimit=1):
    if (protID[0] == "P"):
        protID = protID[1:]
    infos = proteingraph(protID, limit, chemlimit)
    json_file = infos[0]
    protInteractions = listProteinInteractions(protID)
    chemInteractions = listProteinChemicalInteractions(protID)
    allList = []
    allList.extend(protInteractions)
    allList.extend(chemInteractions)
    allInteractions = sorted(allList, key=lambda k: k['combined_score'], reverse=True)
    prot = Proteins.objects.get(protein_id=protID).annotation
    if prot.__contains__(";"):
        name, rest = prot.split(";", 1)
    else:
        name = prot

    template = "stringdb/index.html"
    return render_queryset_to_response(
        request=request,
        template=template,
        data={
            'protid': protID,
            'prot_name': Proteins.objects.get(protein_id=protID).preferred_name,
            'protInteracts': json.dumps(protInteractions),
            'chemInteracts': json.dumps(chemInteractions),
            'allInteracts': json.dumps(allInteractions),
            'name': name,
            'list': ["Homology", "Experiment", "Database", "Textmining",
                     "Genfusion", "Coocurence", "Neighborhood", "Coexpression"],
            'json': json_file,
            'limit': infos[1],
            'protlimit': limit,
            'chemlimit': chemlimit
        }
    )


def index(request):
    proteins = Proteins.objects.filter(species_id=1148)
    print "Schaue nach"
    linkcount_max = 0
    linkcount_min = 1000000
    for protein in proteins:
        counter = NodeNodeLinks.objects.filter(node_id_a=protein.protein_id).count()
        if counter > linkcount_max:
            linkcount_max = counter
        elif counter < linkcount_min:
            linkcount_min = counter

    proteincount = Proteins.objects.count()
    count = proteincount
    return HttpResponse("Min: " + str(linkcount_min) + " Max: " + str(linkcount_max))


''' Getting Interaction Information of an request in an actual showed Network '''


def checkRequest(request):
    data = {}
    selectedID = request.GET["requestID"]
    if selectedID[0] == "P":
        data = checkProteinInteractionRequest(selectedID[1:])
    elif selectedID[0] == "C":
        data = checkChemicalInteractionRequest(selectedID)
    else:
        data = checkProteinInteractionRequest(selectedID)
    template = "stringdb/index.html"
    return render_queryset_to_response(
        request=request,
        template=template,
        data=data
    )


def checkProteinInteractionRequest(protid):
    infos = proteingraph(protid, 10, 1)
    json_file = infos[0]
    protInteractions = listProteinInteractions(protid)
    print "ChemInfos"
    chemInteractions = listProteinChemicalInteractions(protid)
    allList = []
    allList.extend(protInteractions)
    allList.extend(chemInteractions)
    allInteractions = sorted(allList, key=lambda k: k['combined_score'], reverse=True)

    prot = Proteins.objects.get(protein_id=protid).annotation
    if prot.__contains__(";"):
        name, rest = prot.split(";", 1)
    else:
        name = prot
    data = {
        'protid': "P" + str(protid),
        'prot_name': Proteins.objects.get(protein_id=protid).preferred_name,
        'protInteracts': json.dumps(protInteractions).replace("'", "\'"),
        'chemInteracts': json.dumps(chemInteractions).replace("'", "\'"),
        'allInteracts': json.dumps(allInteractions).replace("'", "\'"),
        'name': name,
        'list': ["Homology", "Experiment", "Database", "Textmining",
                 "Genfusion", "Coocurence", "Neighborhood", "Coexpression"],
        'json': json_file,
        'limit': infos[1],
        'protlimit': 10,
        'chemlimit': 1
    }
    return data


def checkChemicalInteractionRequest(chemid):
    infos = chemgraph(chemid, 10, 1)
    json_file = infos[0]
    protInteractions = listChemicalProteinInteractions(chemid)
    print "got protlist"
    chemInteractions = listChemicalChemicalInteractions(chemid)
    print "got chemlist"
    allList = []
    allList.extend(protInteractions)
    allList.extend(chemInteractions)
    allInteractions = sorted(allList, key=lambda k: k['combined_score'], reverse=True)

    chem = Chemicals.objects.get(chemical=chemid).name
    if chem.__contains__(";"):
        name, rest = chem.split(";", 1)
    else:
        name = chem
    data = {
        'protid': chemid,
        'prot_name': chemid,
        'protInteracts': json.dumps(protInteractions),
        'chemInteracts': json.dumps(chemInteractions),
        'allInteracts': json.dumps(allInteractions),
        'name': name,
        'list': ["Homology", "Experiment", "Database", "Textmining",
                 "Genfusion", "Coocurence", "Neighborhood", "Coexpression"],
        'json': json_file,
        'limit': infos[1],
        'protlimit': 10,
        'chemlimit': 1
    }
    return data


''' Getting only the JSON-File for updating the svg Element '''


@ajax_required
def onlyProtGraph(request):
    infos = proteingraph(request.GET["protid"], request.GET["amount"], request.GET["chemAmount"])
    return HttpResponse(infos[0])


''' Loading the information for the given ProteinID with the given amount of neighbours and chemicals'''


def proteingraph(protein_id, limit, chemlimit):
    if protein_id[0] == "P":
        protein_id = protein_id[1:]
    maxhood = 0
    minhood = 10000
    maxScore = 0
    minScore = 10000
    # print protein_id
    protein = Proteins.objects.get(protein_id=protein_id)
    G = nx.Graph()
    interactions = NodeNodeLinks.objects.filter(node_id_a=protein.protein_id).order_by('-combined_score')[:limit]
    # interactions = NodeNodeLinks.objects.filter(node_id_a=prot.protein_id)
    G.add_node("P" + str(protein.protein_id),
               name=protein.preferred_name,
               protein=1,
               hood=0,
               geneid=protein.protein_external_id,
               clicked=1,
               checkbox=1,
               isselected=0,
               annotation=re.sub("'", "", protein.annotation),
               selectvis=1)
    for i in interactions:
        linked = i.node_id_b
        linkedprot = Proteins.objects.get(protein_id=linked)
        hoodsize = len(NodeNodeLinks.objects.filter(node_id_a=linked))
        G.add_node("P" + str(linked),
                   name=linkedprot.preferred_name,
                   protein=1,
                   hood=hoodsize,
                   geneid=linkedprot.protein_external_id,
                   clicked=1,
                   checkbox=1,
                   isselected=0,
                   annotation=re.sub("'", "", linkedprot.annotation),
                   selectvis=1)
        scores = getInteractType(
            NodeNodeLinks.objects.get(node_id_a=protein.protein_id, node_id_b=linked).evidence_scores)
        G.add_edge("P" + str(protein.protein_id),
                   "P" + str(linked),
                   score=i.combined_score,
                   Homology=scores["Homology"],
                   Experiment=scores["Experiment"],
                   Database=scores["Database"],
                   Textmining=scores["Textmining"],
                   Genfusion=scores["Genfusion"],
                   Coocurence=scores["Coocurence"],
                   Neighborhood=scores["Neighborhood"],
                   Coexpression=scores["Coexpression"],
                   clicked=1,
                   checkbox=1,
                   protein=1,
                   selectvis=1)
        maxhood = hoodsize if maxhood < hoodsize else maxhood
        minhood = hoodsize if minhood > hoodsize else minhood
        maxScore = i.combined_score if maxScore < i.combined_score else maxScore
        minScore = i.combined_score if minScore > i.combined_score else minScore

    ''' Looking at each Protein and connecting them, if they interact'''
    for i in interactions:
        for a in interactions:
            if not a.node_id_b == i.node_id_b:
                if not (G.has_edge("P" + str(a.node_id_b), "P" + str(i.node_id_b))):
                    try:
                        link = NodeNodeLinks.objects.get(node_id_a=i.node_id_b, node_id_b=a.node_id_b)
                        scores = getInteractType(link.evidence_scores)
                        G.add_edge("P" + str(i.node_id_b),
                                   "P" + str(a.node_id_b),
                                   score=link.combined_score,
                                   Homology=scores["Homology"],
                                   Experiment=scores["Experiment"],
                                   Database=scores["Database"],
                                   Textmining=scores["Textmining"],
                                   Genfusion=scores["Genfusion"],
                                   Coocurence=scores["Coocurence"],
                                   Neighborhood=scores["Neighborhood"],
                                   Coexpression=scores["Coexpression"],
                                   clicked=1,
                                   checkbox=1,
                                   protein=1,
                                   selectvis=1)
                        maxScore = link.combined_score if maxScore < link.combined_score else maxScore
                        minScore = link.combined_score if minScore > link.combined_score else minScore
                    except ObjectDoesNotExist:
                        pass

    data = findChems(G, chemlimit, maxScore, minScore, minhood, maxhood)
    D = data["graph"]
    maxScore_new = data["maxScore"]
    minScore_new = data["minScore"]
    maxhood_new = data["maxhood"]
    minhood_new = data["minhood"]
    # d = json_graph.node_link_data(G)
    d = json_graph.node_link_data(D)
    json_file = json.dumps(d)
    json_file = calcColourRange(json_file, maxhood_new, minhood_new, maxScore_new, minScore_new)
    data = [json_file, D.number_of_nodes()]
    return data


''' Loading the information for the given ChemicalID with the given amount of neighbours and chemicals'''

# TODO Accelerate Search
def chemgraph(chem_id, limit, chemlimit):
    maxhood = 0
    minhood = 10000
    maxScore = 0
    minScore = 10000
    G = nx.Graph()
    print "looking at chemical"
    chemical = Chemicals.objects.get(chemical=chem_id)
    interactions = ProteinChemicalLinks.objects.filter(chemical=chem_id, protein__regex="1148\.*").order_by('-combined_score')[:limit]
    print "got interacts"
    interactions2 = []
    calledProts = []
    G.add_node(chem_id,
               name=chemical.name,
               protein=0,
               hood=0,
               clicked=1,
               checkbox=1,
               isselected=0,
               selectvis=0)
    for i in interactions:
        linkedprot = Proteins.objects.get(protein_external_id=i.protein)
        if not linkedprot.preferred_name in calledProts:
            calledProts.append(linkedprot.preferred_name)
            interactions2.append(i)
            hoodsize = len(NodeNodeLinks.objects.filter(node_id_a=linkedprot.protein_id))
            if not G.has_node("P" + str(linkedprot.protein_id)):
                G.add_node("P" + str(linkedprot.protein_id),
                           name=linkedprot.preferred_name,
                           protein=1,
                           hood=hoodsize,
                           geneid=linkedprot.protein_external_id,
                           clicked=1,
                           checkbox=1,
                           isselected=0,
                           annotation=re.sub("'", "", linkedprot.annotation),
                           selectvis=1)

                G.add_edge(chem_id, "P" + str(linkedprot.protein_id),
                           score=i.combined_score,
                           clicked=1,
                           checkbox=1,
                           selectvis=0)
                maxhood = hoodsize if maxhood < hoodsize else maxhood
                minhood = hoodsize if minhood > hoodsize else minhood
                maxScore = i.combined_score if maxScore < i.combined_score else maxScore
                minScore = i.combined_score if minScore > i.combined_score else minScore
    print "looked at chem interaction"
    ''' Looking at each Protein and connecting them, if they interact'''
    for i in interactions2:
        prot_i = Proteins.objects.get(protein_external_id=i.protein)
        for a in interactions2:
            prot_a = Proteins.objects.get(protein_external_id=a.protein)
            if not prot_a.protein_id == prot_i.protein_id:
                if not (G.has_edge("P" + str(prot_a.protein_id), "P" + str(prot_i.protein_id))):
                    try:
                        link = NodeNodeLinks.objects.get(node_id_a=prot_i.protein_id, node_id_b=prot_a.protein_id)
                        scores = getInteractType(link.evidence_scores)
                        G.add_edge("P" + str(prot_i.protein_id), "P" + str(prot_a.protein_id),
                                   score=link.combined_score,
                                   Homology=scores["Homology"],
                                   Experiment=scores["Experiment"],
                                   Database=scores["Database"],
                                   Textmining=scores["Textmining"],
                                   Genfusion=scores["Genfusion"],
                                   Coocurence=scores["Coocurence"],
                                   Neighborhood=scores["Neighborhood"],
                                   Coexpression=scores["Coexpression"],
                                   clicked=1,
                                   checkbox=1,
                                   selectvis=1,
                                   protein=1)
                        maxScore = link.combined_score if maxScore < link.combined_score else maxScore
                        minScore = link.combined_score if minScore > link.combined_score else minScore
                    except ObjectDoesNotExist:
                        pass
    print "looked at proteins"

    data = findChems(G, chemlimit, maxScore, minScore, minhood, maxhood)
    print "found chems for proteins"
    D = data["graph"]
    #d = json_graph.node_link_data(G)
    d = json_graph.node_link_data(D)
   #maxScore_new = maxScore
   #minScore_new = minScore
   #maxhood_new = maxhood
   #minhood_new = minhood
    maxScore_new = data["maxScore"]
    minScore_new = data["minScore"]
    maxhood_new = data["maxhood"]
    minhood_new = data["minhood"]
    json_file = json.dumps(d)
    json_file = calcColourRange(json_file, maxhood_new, minhood_new, maxScore_new, minScore_new)
    print "calc color", G.number_of_nodes()

    data = [json_file, G.number_of_nodes()]
    return data


''' Showing the selected Chemical ID with the given amount of Protein and Chemical neighbours from a AJAX Request'''


@ajax_required
def onlyChemGraph(request):
    json = chemgraph(request.GET["chemid"], request.GET["amount"], request.GET["chemAmount"])
    return HttpResponse(json[0])


''' Calculation of the ColourRange for the drawn network and legend and adding it to the JSON-File'''


def calcColourRange(jsonFile, maxHood, minHood, maxScore, minScore):
    colorrangeNode = (maxHood - minHood) / 4
    colorrangeLink = (maxScore - minScore) / 4

    linkcolor = '"linkcolor": [{"color": "red", "value": ' + str(minScore) + '}, ' \
                                                                             '{"color": "yellow", "value": ' + str(
        minScore + colorrangeLink) + '}, ' \
                                     '{"color": "green", "value": ' + str(maxScore - colorrangeLink) + '}, ' \
                                                                                                       '{"color": "navy", "value": ' + str(
        maxScore) + '}]'
    nodecolor = '"nodecolor": [{"color": "black", "value": 0}, ' \
                '{"color": "blue", "value": ' + str(minHood) + '}, ' \
                                                               '{"color": "green", "value": ' + str(
        minHood + colorrangeNode) + '}, ' \
                                    '{"color": "yellow", "value": ' + str(maxHood - colorrangeNode) + '}, ' \
                                                                                                      '{"color": "red", "value": ' + str(
        maxHood) + '}]'
    searchString = '"Experiment":'
    jsonFile = jsonFile.replace(searchString, '"filterscore": 0, ' + searchString)
    jsonFile = jsonFile[:-1] + ', ' + linkcolor + ', ' + nodecolor + '}'
    return jsonFile


''' Getting a List of all Interacts of given ProteinID '''


def listProteinInteractions(protID):
    interactionList = []
    interactions = NodeNodeLinks.objects.filter(node_id_a=protID).order_by("-combined_score")
    for interact in interactions:
        prot = Proteins.objects.get(protein_id=interact.node_id_b)
        protInfo = {}
        name = "P" + str(prot.protein_id)
        protInfo["interactname"] = name
        protInfo["evidence"] = getInteractType(interact.evidence_scores)
        protInfo["preferred_name"] = prot.preferred_name
        protInfo["annotation"] = prot.annotation
        protInfo["combined_score"] = interact.combined_score

        if protInfo["annotation"].__contains__(";"):
            protInfo["annotation"], rest = protInfo["annotation"].split(";", 1)

        interactionList.append(protInfo)
    return interactionList


# TODO changeList Items
def listChemicalChemicalInteractions(chemID):
    interactionList = []
    interactions = ChemicalChemicalLinksDetailed.objects.filter(chemical1=chemID).order_by(
        "-combined_score")
    print "interactions", len(interactions), chemID
    listchems = []
    for interact in interactions:
        try:
            chem = Chemicals.objects.get(chemical=interact.chemical2)
            chemInfo = {}
            chemInfo["interactname"] = chem.chemical
            chemInfo["evidence"] = getChemChemSource(interact)
            chemInfo["preferred_name"] = chem.name
            chemInfo["annotation"] = chemInfo["preferred_name"]
            #chemInfo["formula"] = chem.smiles_string
            chemInfo["combined_score"] = interact.combined_score
            interactionList.append(chemInfo)
        except ObjectDoesNotExist:
            chemInfo = {}
    return interactionList


# ToDo change Listitems
def listProteinChemicalInteractions(protID):
    interactionList = []
    protein = Proteins.objects.get(protein_id=protID).protein_external_id
    interactions = ProteinChemicalLinksDetailed.objects.filter(protein=protein).order_by("-combined_score")
    print "Interactions", len(interactions)
    for interact in interactions:
        chemInfo = {}
        try:
            chem = Chemicals.objects.get(chemical=interact.chemical)

            chemInfo["interactname"] = chem.chemical
            chemInfo["evidence"] = getProteinChemSource(interact.protein, interact.chemical)
            chemInfo["preferred_name"] = chem.name
            chemInfo["annotation"] = chemInfo["preferred_name"]
            #chemInfo["formula"] = chem.smiles_string
            chemInfo["combined_score"] = interact.combined_score
            interactionList.append(chemInfo)
        except ObjectDoesNotExist:
            chemInfo = {}

    return interactionList


def listChemicalProteinInteractions(chemID):
    interactionList = []
    print "looking for interacts"
    interactions = ProteinChemicalLinksDetailed.objects.filter(chemical=chemID, protein__regex="1148\.*").order_by("-combined_score")
    print "got interacts", len(interactions)
    prots = Proteins.objects.filter(protein_external_id__in=map(lambda x: x.protein, interactions))
    for interact, prot in zip(interactions, prots):
        #prot = Proteins.objects.get(protein_external_id=interact.protein)
        protInfo = {}
        name = "P" + str(prot.protein_id)
        protInfo["interactname"] = name
        scoreData = {
            "Experiment": interact.experimental,
            "Database": interact.database,
            "Textmining": interact.textmining,
        }
        protInfo["evidence"] = scoreData
        protInfo["preferred_name"] = prot.preferred_name
        protInfo["annotation"] = prot.annotation
        protInfo["combined_score"] = interact.combined_score

        if protInfo["annotation"].__contains__(";"):
            protInfo["annotation"], rest = protInfo["annotation"].split(";", 1)

        interactionList.append(protInfo)
    print "return list"
    return interactionList


''' Getting the single calculated Scores for out of the StringDB'''


def getInteractType(scores):
    interactType = {"Homology": 0,
                    "Experiment": 0,
                    "Database": 0,
                    "Textmining": 0,
                    "Genfusion": 0,
                    "Coocurence": 0,
                    "Neighborhood": 0,
                    "Coexpression": 0,
                    "CoexpressionA": 0,
                    "CoexpressionB": 0,
                    "NeighborhoodA": 0,
                    "NeighborhoodB": 0,
                    "NeighborhoodC": 0,
                    "CoocurenceA": 0,
                    "CoocurenceB": 0,
                    "GenfusionA": 0,
                    "GenfusionB": 0,
                    "TextminingA": 0,
                    "TextminingB": 0,
                    "DatabaseA": 0,
                    "DatabaseB": 0,
                    "ExperimentA": 0,
                    "ExperimentB": 0,
    }

    for score in scores:
        scoretype = int(score[0])
        if scoretype == 5:
            interactType["Homology"] = score[1]
        elif scoretype == 6:
            interactType["CoexpressionA"] = score[1]
        elif scoretype == 7:
            interactType["CoexpressionB"] = score[1]
        elif scoretype == 8:
            interactType["ExperimentA"] = score[1]
        elif scoretype == 9:
            interactType["ExperimentB"] = score[1]
        elif scoretype == 10:
            interactType["DatabaseA"] = score[1]
        elif scoretype == 11:
            interactType["DatabaseB"] = score[1]
        elif scoretype == 12:
            interactType["TextminingA"] = score[1]
        elif scoretype == 13:
            interactType["TextminingB"] = score[1]
        elif scoretype == 1:
            interactType["NeighborhoodA"] = score[1]
        elif scoretype == 2:
            interactType["NeighborhoodB"] = score[1]
        elif scoretype == 14:
            interactType["NeighborhoodC"] = score[1]
        elif scoretype == 3:
            interactType["GenfusionA"] = score[1]
        elif scoretype == 15:
            interactType["GenfusionB"] = score[1]
        elif scoretype == 4:
            interactType["CoocurenceA"] = score[1]
        elif scoretype == 16:
            interactType["CoocurenceB"] = score[1]

    interactType["Coexpression"] = (1 - (
        (1 - interactType["CoexpressionA"] / 1000.0) * (1 - interactType["CoexpressionB"] / 1000.0))) * 1000.0
    interactType["Experiment"] = (1 - (
        (1 - interactType["ExperimentA"] / 1000.0) * (1 - interactType["ExperimentB"] / 1000.0))) * 1000.0
    interactType["Database"] = (1 - (
        (1 - interactType["DatabaseA"] / 1000.0) * (1 - interactType["DatabaseB"] / 1000.0))) * 1000.0
    interactType["Textmining"] = (1 - (
        (1 - interactType["TextminingA"] / 1000.0) * (1 - interactType["TextminingB"] / 1000.0))) * 1000.0
    interactType["Neighborhood"] = (1 - (
        (1 - interactType["NeighborhoodA"] / 1000.0) * (1 - interactType["NeighborhoodB"] / 1000.0) * (
            1 - interactType["NeighborhoodC"] / 1000.0))) * 1000.0
    interactType["Genfusion"] = (1 - (
        (1 - interactType["GenfusionA"] / 1000.0) * (1 - interactType["GenfusionB"] / 1000.0))) * 1000.0
    interactType["Coocurence"] = (1 - (
        (1 - interactType["CoocurenceA"] / 1000.0) * (1 - interactType["CoocurenceB"] / 1000.0))) * 1000.0

    itemsToRemove =[
        "CoexpressionA",
        "CoexpressionB",
        "NeighborhoodA",
        "NeighborhoodB",
        "NeighborhoodC",
        "CoocurenceA",
        "CoocurenceB",
        "GenfusionA",
        "GenfusionB",
        "TextminingA",
        "TextminingB",
        "DatabaseA",
        "DatabaseB",
        "ExperimentA",
        "ExperimentB"
    ]
    for item in itemsToRemove:
        del interactType[item]
    return interactType


''' Finding Chemicals for each Protein in showed Network.
    Finding min and max Values are deactivated -> Problems during Visualization'''


def findChems(graph, limit, maxScore, minScore, minhood, maxhood):
    graph2 = graph.copy()
    for node in graph2:
        if graph.node[node]["protein"] == 1:
            chemInteracts = ProteinChemicalLinks.objects.filter(protein=graph.node[node]["geneid"]).order_by(
                '-combined_score')[:limit]
            for chemInteract in chemInteracts:
                name = ""
                try:
                    name = Chemicals.objects.get(chemical=chemInteract.chemical).name
                except ObjectDoesNotExist:
                    name = ""
                if not name == "":
                    if not graph.has_node(chemInteract.chemical):
                        hood = len(ProteinChemicalLinks.objects.filter(chemical=chemInteract.chemical,
                                                                       protein__regex="1148\.*"))
                        # if minhood > hood:
                        #    minhood = hood
                        #if maxhood < hood:
                        #    maxhood = hood
                        graph.add_node(chemInteract.chemical,
                                       name=name,
                                       protein=0,
                                       selectvis=0,
                                       clicked=1,
                                       checkbox=1,
                                       hood=hood)
                    score = chemInteract.combined_score
                    # if score > maxScore:
                    #    maxScore = score
                    #if score < minScore:
                    #    minScore = score
                    ProtChemScoreData = getProteinChemSource(node[1:], chemInteract)
                    graph.add_edge(node,
                                   chemInteract.chemical,
                                   score=score,
                                   clicked=1,
                                   checkbox=1,
                                   Homology=0,
                                   Genfusion=0,
                                   Coocurence=0,
                                   Experiment=ProtChemScoreData["Experiment"],
                                   Database=ProtChemScoreData["Database"],
                                   Textmining=ProtChemScoreData["Textmining"],
                                   protein=0,
                                   selectvis=0)
            for chemInteract1 in chemInteracts:
                for chemInteract2 in chemInteracts:
                    if not graph.has_edge(chemInteract1, chemInteract2):
                        try:
                            cheminteraction = ChemicalChemicalLinksDetailed.objects.get(chemical1=chemInteract1,
                                                                                        chemical2=chemInteract2)
                            chemScoreData = getChemChemSource(cheminteraction)
                            graph.add_edge(chemInteract1.chemical,
                                           chemInteract2.chemical,
                                           Genfusion=0,
                                           Coocurence=0,
                                           Homology=chemScoreData["Homology"],
                                           Experiment=chemScoreData["Experiment"],
                                           Database=chemScoreData["Database"],
                                           Textmining=chemScoreData["Textmining"],
                                           score=cheminteraction.combined_score,
                                           clicked=1,
                                           checkbox=1,
                                           protein=0,
                                           selectvis=0)
                        except ObjectDoesNotExist:
                            pass
    data = {'graph': graph, 'minScore': minScore, 'maxScore': maxScore, 'minhood': minhood, 'maxhood': maxhood}
    return data


def getProteinChemSource(protein, chemical):
    data = {}
    try:
        scores = ProteinChemicalLinksDetailed.objects.get(protein=protein, chemical=chemical)
        data = {
            "Experiment": scores.experimental,
            "Database": scores.database,
            "Textmining": scores.textmining,
        }
    except ObjectDoesNotExist:
        data = {
            "Experiment": 0,
            "Database": 0,
            "Textmining": 0,
        }
    return data


''' Resolve the Information of Chemical Chemical Interaction into Dictonary'''


def getChemChemSource(interaction):
    data = {
        "Experiment": interaction.experimental,
        "Database": interaction.database,
        "Textmining": interaction.textmining,
        "Homology": interaction.similarity
    }
    return data


''' Determining Distance between two given Proteins '''


@ajax_required
def getpath(request):
    jsonfile = open("cyanofactory/static/stringdb/cyano.json", "r")
    jsondata = json.load(jsonfile)
    g = json_graph.node_link_graph(jsondata)
    nProtA = 0
    nProtB = 0
    for node in g.nodes():
        if g.node[node]["name"] == request.GET["protA"]:
            nProtA = node
            break

    for node in g.nodes():
        if g.node[node]["name"] == request.GET["protB"]:
            nProtB = node
            break

    found = True

    try:
        result = nx.all_shortest_paths(g, nProtA, nProtB)
        h = nx.Graph()

        nodeattrList = ["name", "hood", "geneid", "clicked", "checkbox", "isselected", "annotation"]
        edgeattrList = ["score", "Homology", "Experiment", "Database", "Textmining", "Genfusion", "Coocurence",
                        "Neighborhood", "Coexpression", "clicked", "checkbox"]

        for r in result:
            n = g.subgraph(r)
            nattrlist = []
            for attr in nodeattrList:
                nattrlist.append(nx.get_node_attributes(n, attr))
            h.add_nodes_from(n.nodes())
            for nameattr, attr in zip(nodeattrList, nattrlist):
                nx.set_node_attributes(h, nameattr, attr)
            eattrList = []
            for attr in edgeattrList:
                eattrList.append(nx.get_edge_attributes(n, attr))
            h.add_edges_from(n.edges())
            for nameattr, attr in zip(edgeattrList, eattrList):
                nx.set_edge_attributes(h, nameattr, attr)

        d = json_graph.node_link_data(h)
        json_file = json.dumps(d)

    except:
        found = False
    return render_queryset_to_response(
        request=request,
        template="stringdb/search.html",
        data={
            "found": found,
            "json": json_file
        }
    )


def startpathway(request):
    return render_queryset_to_response(
        request=request,
        template="stringdb/search.html",
        data={}
    )


'''
Generating JSON-File for each Protein of the organism, retrieving the 10 best scored interaction partner
Until now: looking just at the combined score
To do: Divide between source of information --> Multigraph
'''


def proteingraph2json(request):
    proteins = Proteins.objects.filter(species_id=1148)
    print len(proteins)
    counter = 0
    for prot in proteins:
        counter += 1
        if not os.path.isfile("cyanofactory/stringdb/protein_network/" + str(prot.protein_id) + ".json"):
            print str(counter) + " - " + str(len(proteins))
            G = nx.Graph()
            interactions = NodeNodeLinks.objects.filter(node_id_a=prot.protein_id).order_by('-combined_score')[:10]
            G.add_node(prot.protein_id, name=prot.preferred_name, hood=0, selectvis=1)
            for i in interactions:
                linked = i.node_id_b
                G.add_node(linked, name=Proteins.objects.get(protein_id=linked).preferred_name, hood=len(
                    NodeNodeLinks.objects.filter(node_id_a=linked).order_by('-combined_score')), selectvis=1)
                G.add_edge(prot.protein_id, linked, score=i.combined_score, opacity=1, selectvis=1)

            ''' Looking at each Protein and connecting them, if they interact'''
            for i in interactions:
                for a in interactions:
                    if not a.node_id_b == i.node_id_b:
                        if not (G.has_edge(a.node_id_b, i.node_id_b)):
                            if NodeNodeLinks.objects.filter(node_id_a=i.node_id_b, node_id_b=a.node_id_b).count() != 0:
                                link = NodeNodeLinks.objects.get(node_id_a=i.node_id_b, node_id_b=a.node_id_b)
                                G.add_edge(i.node_id_b,
                                           a.node_id_b,
                                           score=link.combined_score,
                                           opacity=0.2,
                                           selectvis=1)
            print "Schreibe Datei"
            d = json_graph.node_link_data(G)
            data = {"url": 'cyanofactory/stringdb/protein_network/' + str(prot.protein_id) + '.json'}
            json.dump(d, open(data["url"], 'w'))
            a = prot.preferred_name
        else:
            print prot.preferred_name + " already exists"
    a = len(proteins)
    return HttpResponse("Finished - " + str(a))


'''
Comparing 2 Proteins and checking if they have equal interactions partners
Until now 2 Proteins are related if the half of there interactions partners are equal
'''

# TODO Working on a weight function to include Score and amount of interaction partner per Protein
def compareInteraction(request):
    proteins = Proteins.objects.filter(species_id=1148)
    for prot_a in proteins:
        print str(prot_a.protein_id)
        linked_a = NodeNodeLinks.objects.filter(node_id_a=prot_a.protein_id).order_by('-combined_score')
        size_a = len(linked_a)
        print "Protein A " + str(size_a)
        for prot_b in proteins:
            counter = 0
            if not (prot_a.protein_id == prot_b.protein_id):
                linked_b = NodeNodeLinks.objects.filter(node_id_a=prot_b.protein_id).order_by(
                    '-combined_score')
                size_b = len(linked_b)
                for linka in linked_a:
                    linka_id = linka.node_id_b
                    for linkb in linked_b:
                        linkb_id = linkb.node_id_b
                        if linka_id == linkb_id:
                            counter += 1
                            break
                if size_b != 0 and size_a != 0:
                    if (float(counter) / float(size_a)) >= (float(counter) / float(size_b)) and counter != 0:
                        print str(counter) + " " + str(size_a) + " " + str(size_b)
                        print "MATCH"
                        equal = (float(counter) / float(size_a)) / (float(counter) / float(size_b))
                        print str(equal)
                        pc = ProteinComparison
                        pc.objects.get_or_create(protein_a=prot_a.protein_id, protein_b=prot_b.protein_id,
                                                 equal_score=equal)
    return HttpResponse("Done")


'''
Getting the Whole protein_network of the organism.
Needs some time and space
'''


def getwholenetwork():
    print "Running"
    proteins = Proteins.objects.filter(species_id=1148)
    actual = 0
    for protein in proteins:
        if not os.path.isfile("json/" + str(protein.protein_id) + ".json"):
            G = nx.Graph()
            actual += 1
            interactions = NodeNodeLinks.objects.filter(node_id_a=protein.protein_id).order_by('-combined_score')
            if not (G.has_node(protein.protein_id)):
                G.add_node(protein.protein_id,
                           name=protein.preferred_name,
                           hood=0,
                           geneid=protein.protein_external_id,
                           clicked=1,
                           checkbox=1,
                           isselected=0,
                           annotation=re.sub("'", "", protein.annotation),
                           selectvis=1)
            for i in interactions:
                linked = i.node_id_b
                linkedprot = Proteins.objects.get(protein_id=linked)
                hoodsize = len(NodeNodeLinks.objects.filter(node_id_a=linked))
                if not (G.has_node(linked)):
                    G.add_node(linked, name=linkedprot.preferred_name, hood=hoodsize,
                               geneid=linkedprot.protein_external_id, clicked=1, checkbox=1,
                               isselected=0, annotation=re.sub("'", "", protein.annotation), selectvis=1)
                scores = getInteractType(
                    NodeNodeLinks.objects.get(node_id_a=protein.protein_id, node_id_b=linked).evidence_scores)
                if not (G.has_edge(linked, protein.protein_id)):
                    G.add_edge(protein.protein_id,
                               linked,
                               score=i.combined_score,
                               Homology=scores["Homology"],
                               Experiment=scores["Experiment"],
                               Database=scores["Database"],
                               Textmining=scores["Textmining"],
                               Genfusion=scores["Genfusion"],
                               Coocurence=scores["Coocurence"],
                               Neighborhood=scores["Neighborhood"],
                               Coexpression=scores["Coexpression"],
                               clicked=1,
                               checkbox=1,
                               selectvis=1,
                               protein=1)

            d = json_graph.node_link_data(G)
            json.dump(d, open("json/" + str(protein.protein_id) + '.json', 'w'))

    nodeattrList = ["name", "hood", "geneid", "clicked", "checkbox", "isselected", "annotation"]
    edgeattrList = ["score", "Homology", "Experiment", "Database", "Textmining", "Genfusion", "Coocurence",
                    "Neighborhood", "Coexpression", "clicked", "checkbox"]
    g = nx.Graph()
    for afile in os.listdir("json"):
        print afile
        afile = open("json/" + afile, "r")
        adata = json.load(afile)
        ag = json_graph.node_link_graph(adata)

        nattrlist = []
        for attr in nodeattrList:
            nattrlist.append(nx.get_node_attributes(ag, attr))
        g.add_nodes_from(ag.nodes())
        for nameattr, attr in zip(nodeattrList, nattrlist):
            nx.set_node_attributes(g, nameattr, attr)

        eattrList = []
        for attr in edgeattrList:
            eattrList.append(nx.get_edge_attributes(ag, attr))
        g.add_edges_from(ag.edges())
        for nameattr, attr in zip(edgeattrList, eattrList):
            nx.set_edge_attributes(g, nameattr, attr)

        d = json_graph.node_link_data(g)
        output = open("json/cyano.json", "w")
        output.write(json.dumps(d))
    return "finish"