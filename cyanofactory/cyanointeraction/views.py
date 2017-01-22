from django.core.exceptions import ObjectDoesNotExist
from django.http import HttpResponse
from django.shortcuts import render_to_response
from cyano.decorators import ajax_required
from cyanointeraction.models import Proteins, NodeNodeLinks, ProteinsOrthgroups, Chemicals, ProteinChemicalLinks, \
    ChemicalChemicalLinks, ProteinChemicalLinksDetailed, ChemicalChemicalLinksDetailed, StitchNodeNode
import networkx as nx
import json
from networkx.readwrite import json_graph
from cyano.models import ProteinComparison
from cyano.helpers import render_queryset_to_response
import re
import math
import time
import urllib


def checkInteraction(request, protID, protlimit=10, chemlimit=1):
    """
    Getting the interaction information for a selected Protein. Information Contains a Graph of linked Proteins and
    Chemicals.
    @param request:
    @param protID: ID of the selected Protein. It can contain a "P" as first Character or just die ID, which is
    given by the STRING DB
    @param protlimit: maximal amount of Proteins which should be looked at
    @param chemlimit: maximal amount of Chemicals which should be looked at
    @return: render_queryset_to_response which contains the request, the used template for the HTML output and also the
    information of the interaction:
        protid: ProteinID
        prot_name: Proteinname, given by STRING DB
        protInteracts: Protein Interaction Partners
        chemInteracts: Chemical Interaction Partners
        name: preferred_name
        list: List of Interaction possibilities
        json: JSON-File which contains the Interaction-Graph
        limit: number of nodes contained in the Graph
        protlimit: maximal amount of Proteins which should be looked at
        chemlimit: maximal amount of Chemicals which should be looked at
    """
    if protID[0] == "P":
        protID = protID[1:]
    infos = proteingraph(protID, protlimit, chemlimit)
    json_file = infos[0]
    protInteractions = listProteinInteractions(protID)
#    chemInteractions = listProteinChemicalInteractions(protID)
    allList = []
    allList.extend(protInteractions)
#    allList.extend(chemInteractions)
    allInteractions = sorted(allList, key=lambda k: k['combined_score'], reverse=True)
    prot = Proteins.objects.get(protein_id=protID).annotation
    if prot.__contains__(";"):
        name, rest = prot.split(";", 1)
    else:
        name = prot

    template = "cyanointeraction/index.html"
    return render_queryset_to_response(
        request=request,
        template=template,
        data={
            'protid': protID,
            'prot_name': Proteins.objects.get(protein_id=protID).preferred_name,
            'protInteracts': json.dumps(protInteractions),
#            'chemInteracts': json.dumps(chemInteractions),
            'allInteracts': json.dumps(allInteractions),
            'name': name,
            'list': ["Homology", "Experiment", "Database", "Textmining",
                     "Genfusion", "Coocurence", "Neighborhood", "Coexpression"],
            'json': json_file,
            'limit': infos[1],
            'protlimit': protlimit,
            'chemlimit': chemlimit
        }
    )


def index(request):
    proteins = Proteins.objects.filter(species_id=1148)
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


def checkRequest(request):
    """
    Getting the Information for a specified Element from an Interaction.
    It distinguish between a request from a Protein or Chemical and call the appropriate checkInterationRequest.
    @param request: request containing a "requestID"
    @return: render_queryset_to_response with data from checkChemicalInteractionRequest or
    checkProteinInteractionRequest
    """
    data = {}
    selectedID = request.GET.get("requestID", "")
    if len(selectedID) > 0 and selectedID[0] == "P":
        data = checkProteinInteractionRequest(selectedID[1:])
    elif len(selectedID) > 0 and selectedID[0] == "C":
        data = checkChemicalInteractionRequest(int(selectedID[1:]))
    else:
        data = checkProteinInteractionRequest(selectedID)
    template = "cyanointeraction/index.html"
    return render_queryset_to_response(
        request=request,
        template=template,
        data=data
    )


def checkProteinInteractionRequest(protid):
    """
    Generating the data for an Protein-Protein-Interaction. Contains the information of both, Interaction of Proteins
    and Chemicals
    @param protid: Protein ID which is used by the STRING DB
    @return: A Dictionary with all relevant Data:
        protid: ProteinID
        prot_name: Proteinname, given by STRING DB
        protInteracts: Protein Interaction Partners
        chemInteracts: Chemical Interaction Partners
        name: preferred_name
        list: List of Interaction possibilities
        json: JSON-File which contains the Interaction-Graph
        limit: number of nodes contained in the Graph
        protlimit: maximal amount of Proteins which should be looked at
        chemlimit: maximal amount of Chemicals which should be looked at
    """
    infos = proteingraph(protid, 10, 1)
    json_file = infos[0]
    protInteractions = listProteinInteractions(protid)
#    chemInteractions = listProteinChemicalInteractions(protid)
    allList = []
    allList.extend(protInteractions)
#    allList.extend(chemInteractions)
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
#        'chemInteracts': json.dumps(chemInteractions).replace("'", "\'"),
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



def checkSTITCH(chemid):
    """
    Generating the data for an Chemical-Protein-Interaction. Contains the information of both, Interaction of Proteins
    and Chemicals
    @param chemid: CHEMID which is used by the STITCH DB without leading "-"
    @return: A Dictionary with all relevant Data:
        protid: ProteinID
        prot_name: Proteinname, given by STRING DB
        protInteracts: Protein Interaction Partners
        chemInteracts: Chemical Interaction Partners
        name: preferred_name
        list: List of Interaction possibilities
        json: JSON-File which contains the Interaction-Graph
        limit: number of nodes contained in the Graph
        protlimit: maximal amount of Proteins which should be looked at
        chemlimit: maximal amount of Chemicals which should be looked at
    """

    infos = stitchgraph(chemid, 10, 1)
    json_file = infos[0]
    protInteractions = listSTITCHProteinInteractions(chemid)
    chemInteractions = listSTITCHChemicalInteractions(chemid)
    allList = []
    allList.extend(protInteractions)
    allList.extend(chemInteractions)
    allInteractions = sorted(allList, key=lambda k: k['combined_score'], reverse=True)

    chem = Chemicals.objects.get(chemical_id=-chemid).short_name
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

def testSTITCH(a):
    checkSTITCH(171549453)

def checkChemicalInteractionRequest(chemid):
    """
    Generating the data for an Chemical-Protein-Interaction. Contains the information of both, Interaction of Proteins
    and Chemicals
    @param chemid: CHEMID which is used by the STITCH DB
    @return: A Dictionary with all relevant Data:
        protid: ProteinID
        prot_name: Proteinname, given by STRING DB
        protInteracts: Protein Interaction Partners
        chemInteracts: Chemical Interaction Partners
        name: preferred_name
        list: List of Interaction possibilities
        json: JSON-File which contains the Interaction-Graph
        limit: number of nodes contained in the Graph
        protlimit: maximal amount of Proteins which should be looked at
        chemlimit: maximal amount of Chemicals which should be looked at
    """
    infos = chemgraph(chemid, 10, 1)
    json_file = infos[0]
    protInteractions = listChemicalProteinInteractions(chemid)
    chemInteractions = listChemicalChemicalInteractions(chemid)
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


@ajax_required
def onlyProtGraph(request):
    """
    Getting only the Graph of an Protein-Interactions from an ajax request.
    @param request: it must contain a protid, amount and chemAmount
    @return: Graph in JSON-Format
    """
    infos = proteingraph(request.GET["protid"], request.GET["amount"], request.GET["chemAmount"])
    return HttpResponse(infos[0])


def proteingraph(protein_id, protlimit, chemlimit):
    """
    Loading the information for a given Protein ID with a given amount of protein and chemical interaction partner to
    generate an interaction graph.
    @param protein_id: Protein ID given by STRING DB. It can start with a "P" or not
    @param protlimit: maximal amount of protein interaction partners
    @param chemlimit: maximal amount of chemical interaction partners
    @return: List with 2 entries, which contain data from the generated graph:
        First Element - graph in JSON-Format
        Second Element - amount of nodes from the graph
    """
    if protein_id[0] == "P":
        protein_id = protein_id[1:]

    geneList = readGeneList()

    maxhood = 0
    minhood = 10000
    maxScore = 0
    minScore = 10000
    # print protein_id
    protein = Proteins.objects.get(protein_id=protein_id)
    G = nx.Graph()
    interactions = NodeNodeLinks.objects.filter(node_id_a=protein.protein_id).order_by('-combined_score')[:protlimit]
    # interactions = NodeNodeLinks.objects.filter(node_id_a=prot.protein_id)
    changed = 0
    if protein_id in geneList:
        changed = 1
    G.add_node("P" + str(protein.protein_id),
               name=protein.preferred_name,
               protein=1,
               hood=0,
               geneid=protein.protein_id,
               clicked=1,
               checkbox=1,
               isselected=0,
               annotation=re.sub("'", "", protein.annotation),
               selectvis=1)
    for i in interactions:
        linked = i.node_id_b
        if linked in geneList:
            changed = 1
        linkedprot = Proteins.objects.get(protein_id=linked)
        hoodsize = NodeNodeLinks.objects.filter(node_id_a=linked).count()
        G.add_node("P" + str(linked),
                   name=linkedprot.preferred_name,
                   protein=1,
                   hood=hoodsize,
                   geneid=linkedprot.protein_id,
                   clicked=1,
                   checkbox=1,
                   isselected=0,
                   annotation=re.sub("'", "", linkedprot.annotation),
                   selectvis=1)
        scores = getInteractType(i.evidence_scores)
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
                   selectvis=1,
                   snp=changed)
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
                        changedCon=0
                        if a.node_id_b in geneList or i.node_id_a in geneList:
                            changedCon=1
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
                                   selectvis=1,
                                   snp=changedCon)
                        maxScore = link.combined_score if maxScore < link.combined_score else maxScore
                        minScore = link.combined_score if minScore > link.combined_score else minScore
                    except ObjectDoesNotExist:
                        pass
# todo fix Stitch fixing
#    data = findChems(G, chemlimit)
#    data = findSTITCHChems(G, chemlimit)
#    D = data
    d = json_graph.node_link_data(G)
#    d = json_graph.node_link_data(D)
    json_file = json.dumps(d)
    json_file = calcColourRange(json_file, maxhood, minhood, maxScore, minScore)
#    data = [json_file, D.number_of_nodes()]
    data = [json_file, G.number_of_nodes()]
    return data


def stitchgraph(chem_id, protlimit, chemlimit):
    """
    Loading the information for a given Chemical ID with a given amount of protein and chemical interaction partner to
    generate an interaction graph.
    @param chem_id: CHEMICAL ID used in the STITCH DB
    @param protlimit: maximal amount of protein interaction partners
    @param chemlimit: maximal amount of chemical interaction partners
    @return: List with 2 entries, which contain data from the generated graph:
        First Element - graph in JSON-Format
        Second Element - amount of nodes from the graph
    """
    geneList = readGeneList()

    maxhood = 0
    minhood = 10000
    maxScore = 0
    minScore = 10000
    G = nx.Graph()
    chemical = Chemicals.objects.get(chemical_id=-chem_id)
    interactions = StitchNodeNode.objects.filter(node_id_a=-chem_id, node_type_b=1148).order_by('-combined_score')[:protlimit]
    interactions2 = []
    calledProts = []
    G.add_node(chem_id,
               name=chemical.preferred_name,
               protein=0,
               hood=0,
               clicked=1,
               checkbox=1,
               isselected=0,
               selectvis=0)
    for i in interactions:
        linkedprot = Proteins.objects.get(protein_id=i.node_id_b)
        if not linkedprot.preferred_name in calledProts:
            calledProts.append(linkedprot.preferred_name)
            interactions2.append(i)
            change = 0
            if linkedprot.protein_id in geneList:
                change = 1
            hoodsize = NodeNodeLinks.objects.filter(node_id_a=linkedprot.protein_id).count()
            if not G.has_node("P" + str(linkedprot.protein_id)):
                G.add_node("P" + str(linkedprot.protein_id),
                           name=linkedprot.preferred_name,
                           protein=1,
                           hood=hoodsize,
                           geneid=linkedprot.protein_id,
                           clicked=1,
                           checkbox=1,
                           isselected=0,
                           annotation=re.sub("'", "", linkedprot.annotation),
                           selectvis=1)

                G.add_edge(chem_id, "P" + str(linkedprot.protein_id),
                           score=i.combined_score,
                           clicked=1,
                           checkbox=1,
                           selectvis=0,
                           snp=change)
                maxhood = hoodsize if maxhood < hoodsize else maxhood
                minhood = hoodsize if minhood > hoodsize else minhood
                maxScore = i.combined_score if maxScore < i.combined_score else maxScore
                minScore = i.combined_score if minScore > i.combined_score else minScore
    ''' Looking at each Protein and connecting them, if they interact'''
    for i in interactions2:
        prot_i = Proteins.objects.get(protein_id=i.node_id_b)
        for a in interactions2:
            prot_a = Proteins.objects.get(protein_id=a.node_id_b)
            if not prot_a.protein_id == prot_i.protein_id:
                if not (G.has_edge("P" + str(prot_a.protein_id), "P" + str(prot_i.protein_id))):
                    try:
                        link = NodeNodeLinks.objects.get(node_id_a=prot_i.protein_id, node_id_b=prot_a.protein_id)
                        scores = getInteractType(link.evidence_scores)
                        changedCon = 0
                        if prot_i.protein_id in geneList or prot_a.protein_id in geneList:
                            changedCon = 1
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
                                   protein=1,
                                   snp=changedCon)
                        maxScore = link.combined_score if maxScore < link.combined_score else maxScore
                        minScore = link.combined_score if minScore > link.combined_score else minScore
                    except ObjectDoesNotExist:
                        pass

    #data = findSTITCHChems(G, chemlimit)
    #D = data
    d = json_graph.node_link_data(G)
    #d = json_graph.node_link_data(D)
    json_file = json.dumps(d)
    json_file = calcColourRange(json_file, maxhood, minhood, maxScore, minScore)

    data = [json_file, G.number_of_nodes()]
    return data

# TODO Accelerate Search
def chemgraph(chem_id, protlimit, chemlimit):
    """
    Loading the information for a given Chemical ID with a given amount of protein and chemical interaction partner to
    generate an interaction graph.
    @param chem_id: CHEMICAL ID used in the STITCH DB
    @param protlimit: maximal amount of protein interaction partners
    @param chemlimit: maximal amount of chemical interaction partners
    @return: List with 2 entries, which contain data from the generated graph:
        First Element - graph in JSON-Format
        Second Element - amount of nodes from the graph
    """

    geneList = readGeneList()

    maxhood = 0
    minhood = 10000
    maxScore = 0
    minScore = 10000
    G = nx.Graph()
    print(-chem_id)
    chemical = Chemicals.objects.get(chemical_id=-chem_id)
    interactions = ProteinChemicalLinks.objects.filter(chemical=chem_id, protein__startswith="1148").order_by('-combined_score')[:protlimit]
    interactions2 = []
    calledProts = []
    G.add_node(chem_id,
               name=chemical.preferred_name,
               protein=0,
               hood=0,
               clicked=1,
               checkbox=1,
               isselected=0,
               selectvis=0)
    for i in interactions:
        linkedprot = Proteins.objects.get(protein_id=i.protein)
        if not linkedprot.preferred_name in calledProts:
            calledProts.append(linkedprot.preferred_name)
            interactions2.append(i)
            change = 0
            if linkedprot.protein_id in geneList:
                change = 1
            hoodsize = NodeNodeLinks.objects.filter(node_id_a=linkedprot.protein_id).count()
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
                           selectvis=0,
                           snp=change)
                maxhood = hoodsize if maxhood < hoodsize else maxhood
                minhood = hoodsize if minhood > hoodsize else minhood
                maxScore = i.combined_score if maxScore < i.combined_score else maxScore
                minScore = i.combined_score if minScore > i.combined_score else minScore
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
                        changedCon = 0
                        if prot_i.protein_id in geneList or prot_a.protein_id in geneList:
                            changedCon = 1
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
                                   protein=1,
                                   snp=changedCon)
                        maxScore = link.combined_score if maxScore < link.combined_score else maxScore
                        minScore = link.combined_score if minScore > link.combined_score else minScore
                    except ObjectDoesNotExist:
                        pass

    data = findChems(G, chemlimit)
    D = data
    #d = json_graph.node_link_data(G)
    d = json_graph.node_link_data(D)
    json_file = json.dumps(d)
    json_file = calcColourRange(json_file, maxhood, minhood, maxScore, minScore)

    data = [json_file, G.number_of_nodes()]
    return data


@ajax_required
def onlyChemGraph(request):
    """
    Showing the selected Chemical ID with a given amount of Protein and Chemical neighbours from a ajax Request
    @param request: needs the entries chemid, amount and chemAmount
    @return: Graph in JSON-Format
    """
    json_graph = chemgraph(request.GET["chemid"], request.GET["amount"], request.GET["chemAmount"])
    return HttpResponse(json_graph[0])


def calcColourRange(jsonFile, maxHood, minHood, maxScore, minScore):
    """
    Calculation of the ColourRange for the drawn graph and legend and adding it to the JSON-File, which contains the
    graph.
    @param jsonFile: graph in JSON-Format
    @param maxHood: maximal amount of nodes neighbours in graph
    @param minHood: mininmal amount of nodes neighbours in graph
    @param maxScore: maximal interaction Score in graph
    @param minScore: minimal interaction Score in graph
    @return: JSON formatted string containing the graph an the color range for nodes and edges
    """
    colorrangeNode = (maxHood - minHood) / 4
    colorrangeLink = (maxScore - minScore) / 4

    linkcolor = '"linkcolor": [{"color": "red", "value": ' + str(minScore) + '}, ' \
                                '{"color": "yellow", "value": ' + str(minScore + colorrangeLink) + '}, ' \
                                '{"color": "green", "value": ' + str(maxScore - colorrangeLink) + '}, ' \
                                '{"color": "navy", "value": ' + str(maxScore) + '}]'
    nodecolor = '"nodecolor": [{"color": "black", "value": 0}, ' \
                                '{"color": "blue", "value": ' + str(minHood) + '}, ' \
                                '{"color": "green", "value": ' + str(minHood + colorrangeNode) + '}, ' \
                                '{"color": "yellow", "value": ' + str(maxHood - colorrangeNode) + '}, ' \
                                '{"color": "red", "value": ' + str(maxHood) + '}]'
    searchString = '"Experiment":'
    jsonFile = jsonFile.replace(searchString, '"filterscore": 0, ' + searchString)
    jsonFile = jsonFile[:-1] + ', ' + linkcolor + ', ' + nodecolor + '}'
    return jsonFile


def listProteinInteractions(protID):
    """
    Getting all protein interaction partners by a given protein ID and sorted by their combined score
    @param protID: protein ID from STRING DB
    @return: List of dictionaries containing following informations:
        interactname - Protein ID
        evidence = Dictionary of interaction possibilities and there scores
        preferred_name = Name of the Protein
        annotation = additional annotation of a protein
        combined_score = interaction score of protein with selected protein
    """
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

        if ";" in protInfo["annotation"]:
            protInfo["annotation"], rest = protInfo["annotation"].split(";", 1)

        interactionList.append(protInfo)
    return interactionList


def listSTITCHChemicalInteractions(chemID):
    """
    Getting all chemical interaction partners by a given chemical ID and sorted by their combined score
    @param chemID: chemical ID by STITCH DB
    @return: List of dictionaries containing following informations:
        interactname - ID of the chemical
        evidence = Dictionary of interaction possibilities and there scores
        preferred_name = chemical.preferred_name
        annotation = additional annotation of a protein
        combined_score = interaction score of chemical with selected protein
    """
    interactionList = []
    interactions = StitchNodeNode.objects.filter(node_id_a=-chemID, node_id_b__startswith="-").order_by(
        "-combined_score")
    for interact in interactions:
        try:
            chem = Chemicals.objects.get(chemical_id=interact.node_id_b)
            chemInfo = {"interactname": chem.short_name,
                        "evidence": getChemChemSource(interact),
                        "preferred_name": chem.preferred_name,
                        "annotation": "",
                        "combined_score": interact.combined_score
                        }

            #chemInfo["formula"] = chem.smiles_string
            interactionList.append(chemInfo)
        except ObjectDoesNotExist:
            chemInfo = {}
    return interactionList


# TODO changeList Items
def listChemicalChemicalInteractions(chemID):
    """
    Getting all chemical interaction partners by a given chemical ID and sorted by their combined score
    @param chemID: chemical ID by STITCH DB
    @return: List of dictionaries containing following informations:
        interactname - ID of the chemical
        evidence = Dictionary of interaction possibilities and there scores
        preferred_name = chemical.preferred_name
        annotation = additional annotation of a protein
        combined_score = interaction score of chemical with selected protein
    """
    interactionList = []
    interactions = ChemicalChemicalLinksDetailed.objects.filter(chemical1=chemID).order_by(
        "-combined_score")
    for interact in interactions:
        try:
            chem = Chemicals.objects.get(chemical=interact.chemical2)
            chemInfo = {}
            chemInfo["interactname"] = chem.chemical
            chemInfo["evidence"] = getChemChemSource(interact)
            chemInfo["preferred_name"] = chem.name
            chemInfo["annotation"] = ""
            #chemInfo["formula"] = chem.smiles_string
            chemInfo["combined_score"] = interact.combined_score
            interactionList.append(chemInfo)
        except ObjectDoesNotExist:
            chemInfo = {}
    return interactionList


# ToDo change Listitems
def listProteinChemicalInteractions(protID):
    """
    Getting all chemical interaction partners by a given protein ID and sorted by their combined score
    @param protID: protein ID from STRING DB
    @return: List of dictionaries containing following informations:
        interactname - ID of the chemical
        evidence = Dictionary of interaction possibilities and there scores
        preferred_name = chemical.preferred_name
        annotation = additional annotation of a protein
        combined_score = interaction score of chemical with selected protein
    """
    interactionList = []
    listChems = []
    protein = Proteins.objects.get(protein_id=protID).protein_external_id
    interactions = ProteinChemicalLinksDetailed.objects.filter(protein=protein).order_by("-combined_score")
    for interact in interactions:
        chemInfo = {}
        try:
            chem = Chemicals.objects.get(chemical=interact.chemical)
            if not chem.name in listChems:
                chemInfo["interactname"] = chem.chemical
                chemInfo["evidence"] = getProteinChemSource(interact.protein, interact.chemical)
                chemInfo["preferred_name"] = chem.name
                chemInfo["annotation"] = ""
                #chemInfo["formula"] = chem.smiles_string
                chemInfo["combined_score"] = interact.combined_score
                listChems.append(chem.name)
                interactionList.append(chemInfo)
        except ObjectDoesNotExist:
            chemInfo = {}

    return interactionList


def listSTITCHProteinInteractions(chemID):
    """
    Getting all protein interaction partners by a given chemical ID and sorted by their combined score
    @param chemID: chemical ID by STITCH DB
    @return: List of dictionaries containing following informations:
        interactname - Protein ID
        evidence = Dictionary of interaction possibilities and there scores
        preferred_name = Name of the Protein
        annotation = additional annotation of a protein
        combined_score = interaction score of protein with selected chemical
    """
    interactionList = []
    interactions = StitchNodeNode.objects.filter(node_id_a=-chemID, node_type_b=1148).order_by("-combined_score")
    prots = Proteins.objects.filter(protein_id__in=map(lambda x: x.node_id_b, interactions))
    for interact, prot in zip(interactions, prots):
        #prot = Proteins.objects.get(protein_external_id=interact.protein)
        protInfo = {}
        name = "P" + str(prot.protein_id)
        protInfo["interactname"] = name
        scoreData = {}
        """
        scoreData = {
            "Experiment": interact.experimental,
            "Database": interact.database,
            "Textmining": interact.textmining,
        }
        """
        protInfo["evidence"] = scoreData
        protInfo["preferred_name"] = prot.preferred_name
        protInfo["annotation"] = prot.annotation
        protInfo["combined_score"] = interact.combined_score

        if protInfo["annotation"].__contains__(";"):
            protInfo["annotation"], rest = protInfo["annotation"].split(";", 1)

        interactionList.append(protInfo)
    return interactionList


def listChemicalProteinInteractions(chemID):
    """
    Getting all protein interaction partners by a given chemical ID and sorted by their combined score
    @param chemID: chemical ID by STITCH DB
    @return: List of dictionaries containing following informations:
        interactname - Protein ID
        evidence = Dictionary of interaction possibilities and there scores
        preferred_name = Name of the Protein
        annotation = additional annotation of a protein
        combined_score = interaction score of protein with selected chemical
    """
    interactionList = []
    interactions = ProteinChemicalLinksDetailed.objects.filter(chemical=chemID, protein__startswith="1148").order_by("-combined_score")
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
    return interactionList


#ToDO add param
def getInteractType(scores):
    """
    Getting the single calculated Scores out of the StringDB
    @param scores: Object from....
    @return: Dictionary with all score types
    """
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


def findSTITCHChems(graph, chemLimit):
    """
    Finding Chemicals for each Protein in given network.
    !!!
     Finding min and max Values are deactivated -> Problems during Visualization. No max/min calculation for
     chemicals
    !!!
    @param graph: networkx graph
    @param protlimit: maximal amount of Proteins which should be called for each chemical
    @param maxScore: maximal combined_score
    @param minScore: minimal combined_score
    @param minhood: minimal neighborhood of a node
    @param maxhood: maximal neighborhood of a node
    @return: networkx graph with added chemical nodes
    """
    graph2 = graph.copy()
    for node in graph2:
        if graph.node[node]["protein"] == 1:
            chemInteracts = StitchNodeNode.objects.filter(node_id_a=graph.node[node]["geneid"]).order_by(
                '-combined_score')[:chemLimit]
            for chemInteract in chemInteracts:
                name = ""
                try:
                    name = Chemicals.objects.get(chemical_id=chemInteract.node_id_b).name
                except ObjectDoesNotExist:
                    name = ""
                if not name == "":
                    if not graph.has_node(chemInteract.chemical):
                        hood = ProteinChemicalLinks.objects.filter(chemical=chemInteract.chemical,
                                                                       protein__startswith="1148").count()
                        graph.add_node(chemInteract.chemical,
                                       name=name,
                                       protein=0,
                                       selectvis=0,
                                       clicked=1,
                                       checkbox=1,
                                       hood=hood)
                    score = chemInteract.combined_score
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
                            cheminteraction = StitchNodeNode.objects.get(node_id_a=chemInteract1.node_id_a,
                                                                         node_id_b=chemInteract2.node_id_a)
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
    return graph


def findChems(graph, chemLimit):
    """
    Finding Chemicals for each Protein in given network.
    !!!
     Finding min and max Values are deactivated -> Problems during Visualization. No max/min calculation for
     chemicals
    !!!
    @param graph: networkx graph
    @param protlimit: maximal amount of Proteins which should be called for each chemical
    @param maxScore: maximal combined_score
    @param minScore: minimal combined_score
    @param minhood: minimal neighborhood of a node
    @param maxhood: maximal neighborhood of a node
    @return: networkx graph with added chemical nodes
    """
    graph2 = graph.copy()
    for node in graph2:
        if graph.node[node]["protein"] == 1:
            chemInteracts = ProteinChemicalLinks.objects.filter(protein=graph.node[node]["geneid"]).order_by(
                '-combined_score')[:chemLimit]
            for chemInteract in chemInteracts:
                name = ""
                try:
                    name = Chemicals.objects.get(chemical=chemInteract.chemical).name
                except ObjectDoesNotExist:
                    name = ""
                if not name == "":
                    if not graph.has_node(chemInteract.chemical):
                        hood = ProteinChemicalLinks.objects.filter(chemical=chemInteract.chemical,
                                                                       protein__startswith="1148").count()
                        graph.add_node(chemInteract.chemical,
                                       name=name,
                                       protein=0,
                                       selectvis=0,
                                       clicked=1,
                                       checkbox=1,
                                       hood=hood)
                    score = chemInteract.combined_score
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
    return graph


def getProteinChemSource(protein, chemical):
    """
    Collecting single scores of an interaction for an given protein and chemical
    @param protein: external protein ID from STRING/ protein ID from STITCH
    @param chemical: CHEMID from STITCH
    @return: dictionary containing the single scores: Experiment, Database and Textmining
    """
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


def getChemChemSource(interaction):
    """
    Collecting single scores of an interaction for an given chemical chemical interaction
    @param interaction: ChemicalChemicalLinksDetailed object
    @return: dictionary with single interaction scores: Experiment, Database, Textmining, Homology
    """
    data = {
        "Experiment": "",#interaction.experimental,
        "Database": "",#interaction.database,
        "Textmining": "",#interaction.textmining,
        "Homology": "",#interaction.similarity
    }
    return data


@ajax_required
def getpath(request):
    """
    Getting the path between 2 given nodes from an ajax request
    @param request: contains protA and protB. Both are protein ids from STRING
    @return: render_queryset_to_response with defined html page and data dictionary containing an boolean value for
        found of a path an a graph containing the path.
    """
    jsonfile = open("cyanofactory/static/cyanointeraction/cyano.json", "r")
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
        template="cyanointeraction/search.html",
        data={
            "found": found,
            "json": json_file
        }
    )


def startpathway(request):
    return render_queryset_to_response(
        request=request,
        template="cyanointeraction/search.html",
        data={}
    )


'''
Generating JSON-File for each Protein of the organism, retrieving the 10 best scored interaction partner
Until now: looking just at the combined score
To do: Divide between source of information --> Multigraph
'''


def proteingraph2json(request):
    proteins = Proteins.objects.filter(species_id=1148)
    #print len(proteins)
    counter = 0
    for prot in proteins:
        counter += 1
        if not os.path.isfile("cyanofactory/cyanointeraction/protein_network/" + str(prot.protein_id) + ".json"):
            #print str(counter) + " - " + str(len(proteins))
            G = nx.Graph()
            interactions = NodeNodeLinks.objects.filter(node_id_a=prot.protein_id).order_by('-combined_score')[:10]
            G.add_node(prot.protein_id, name=prot.preferred_name, hood=0, selectvis=1)
            for i in interactions:
                linked = i.node_id_b
                G.add_node(linked, name=Proteins.objects.get(protein_id=linked).preferred_name, hood=
                    NodeNodeLinks.objects.filter(node_id_a=linked).order_by('-combined_score').count(), selectvis=1)
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
            print("Schreibe Datei")
            d = json_graph.node_link_data(G)
            data = {"url": 'cyanofactory/cyanointeraction/protein_network/' + str(prot.protein_id) + '.json'}
            json.dump(d, open(data["url"], 'w'))
            a = prot.preferred_name
        else:
            print(prot.preferred_name + " already exists")
    a = proteins.count()
    return HttpResponse("Finished - " + str(a))


'''
Comparing 2 Proteins and checking if they have equal interactions partners
Until now 2 Proteins are related if the half of there interactions partners are equal
'''

# TODO Working on a weight function to include Score and amount of interaction partner per Protein
def compareInteraction(request):
    proteins = Proteins.objects.filter(species_id=1148)
    for prot_a in proteins:
        print(str(prot_a.protein_id))
        linked_a = NodeNodeLinks.objects.filter(node_id_a=prot_a.protein_id).order_by('-combined_score')
        size_a = len(linked_a)
        print("Protein A " + str(size_a))
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
                        print(str(counter) + " " + str(size_a) + " " + str(size_b))
                        print("MATCH")
                        equal = (float(counter) / float(size_a)) / (float(counter) / float(size_b))
                        print(str(equal))
                        pc = ProteinComparison
                        pc.objects.get_or_create(protein_a=prot_a.protein_id, protein_b=prot_b.protein_id,
                                                 equal_score=equal)
    return HttpResponse("Done")


'''
Getting the Whole protein_network of the organism.
Needs some time and space
'''


def getwholenetwork():
    print("Running")
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
        print(afile)
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


def getTimes(ids):
    counter = 0
    procent = 0
    timeList = []
    ids = ids[:len(ids)]
    for id in ids:
        startTime = time.time()
        temp = urllib.urlopen("http://127.0.0.1:8000/cyanointeraction/interaction/" + str(id) + "/")
        temp.read()
        temp.close()
        endTime = time.time()
        workTime = endTime - startTime
        timeList.append(workTime)
        counter += 1
        if counter == len(ids)/10:
            counter = 0
            procent += 10
            print(str(procent)+"% done")
    return timeList


def writeTimeStatistic(timeList, writeFile):
    meantime = 0
    for time in timeList:
        writeFile.write(str(time) + "\n")
        meantime += time

    meantime /= len(timeList)
    writeFile.write("\n" + str(meantime) + "\n")

    std = 0
    for time in timeList:
        std += math.sqrt(math.pow(time - meantime, 2))

    std /= len(timeList) - 1
    writeFile.write("\n" + str(std) + "\n")

def calcTimes(request):
    listProteins = Proteins.objects.filter(species_id=1148)
    listProteinIDs = []
    for prot in listProteins:
        listProteinIDs.append(prot.protein_id)
    listProtTime = getTimes(listProteinIDs)

    outputProtTime = open("protTime2", "w")
    writeTimeStatistic(listProtTime, outputProtTime)
    outputProtTime.close()

   ##listChems = ProteinChemicalLinksDetailed.objects.filter(protein__regex="1148\.*")
   ##listChemIDs = []
   ##for chem in listChems:
   ##    listChemIDs.append(chem.chemical)
   ##listChemTime = getTimes(listChemIDs)

   ##outputChemTime = open("chemTime", "w")
   ##writeTimeStatistic(listChemTime, outputChemTime)
   ##outputChemTime.close()


    return HttpResponse("Done")


import os
from django.conf import settings
def readGeneList():
    # todo check if line ending is correct
    print(settings.STATICFILES_DIRS[0])
    i_handler = open(settings.STATICFILES_DIRS[0]+"/cyanointeraction/genes_with_snps.txt", "r")
    gene_list = [line for line in i_handler]
    return gene_list