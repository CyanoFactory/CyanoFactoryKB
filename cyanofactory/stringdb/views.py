from django.core.exceptions import ObjectDoesNotExist
from django.http import HttpResponse
from django.shortcuts import render_to_response
from cyano.decorators import ajax_required
from stringdb.models import Proteins, NodeNodeLinks, ProteinsOrthgroups
import networkx as nx
import json
from networkx.readwrite import json_graph
import os
from cyano.models import ProteinComparison
from cyano.helpers import render_queryset_to_response
import re


def checkInteraction(request, protID, limit=10):
    json_file = proteingraph(protID, limit)
    interactions = listInteractions(protID)
    prot = Proteins.objects.get(protein_id=protID).annotation
    if prot.__contains__(";"):
        name, rest = prot.split(";", 1)
    else:
        name = prot

    if request.is_ajax():
        template = "stringdb/list_page.html"
    else:
        template = "stringdb/index.html"
    return render_queryset_to_response(
        request=request,
        template=template,
        data={
            'protid': protID,
            'prot_name': Proteins.objects.get(protein_id=protID).preferred_name,
            'interacts': interactions,
            'name': name,
            'list': ["Homology", "Experiment", "Database", "Textmining",
                     "Genfusion", "Coocurence", "Neighborhood", "Coexpression"],
            'json': json_file,
            'limit': limit
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

@ajax_required
def onlygraph(request):
    json = proteingraph(request.GET["protid"], request.GET["amount"])
    return HttpResponse(json)

def proteingraph(protein_id, limit):
    maxhood = 0
    minhood = 10000
    maxScore = 0
    minScore = 1000
    #print protein_id
    protein = Proteins.objects.get(protein_id=protein_id)
    G = nx.Graph()
    interactions = NodeNodeLinks.objects.filter(node_id_a=protein.protein_id).order_by('-combined_score')[:limit]
    #interactions = NodeNodeLinks.objects.filter(node_id_a=prot.protein_id)
    G.add_node(protein.protein_id, name=protein.preferred_name, hood=0, geneid=protein.protein_external_id, clicked=1, checkbox=1, isselected=0)
    for i in interactions:
        linked = i.node_id_b
        hoodsize = len(NodeNodeLinks.objects.filter(node_id_a=linked))
        G.add_node(linked, name=Proteins.objects.get(protein_id=linked).preferred_name, hood=hoodsize, geneid=Proteins.objects.get(protein_id=linked).protein_external_id, clicked=1, checkbox=1, isselected=0)
        scores = getInteractType(
            NodeNodeLinks.objects.get(node_id_a=protein.protein_id, node_id_b=linked).evidence_scores)
        G.add_edge(protein.protein_id, linked, score=i.combined_score,
                   Homology=scores["Homology"], Experiment=scores["Experiment"],
                   Database=scores["Database"], Textmining=scores["Textmining"],
                   Genfusion=scores["Genfusion"], Coocurence=scores["Coocurence"],
                   Neighborhood=scores["Neighborhood"], Coexpression=scores["Coexpression"], clicked=1, checkbox=1)
        maxhood = hoodsize if maxhood < hoodsize else maxhood
        minhood = hoodsize if minhood > hoodsize else minhood
        maxScore = i.combined_score if maxScore < i.combined_score else maxScore
        minScore = i.combined_score if minScore > i.combined_score else minScore

    ''' Looking at each Protein and connecting them, if they interact'''
    for i in interactions:
        for a in interactions:
            if not a.node_id_b == i.node_id_b:
                if not (G.has_edge(a.node_id_b, i.node_id_b)):
                    try:
                        link = NodeNodeLinks.objects.get(node_id_a=i.node_id_b, node_id_b=a.node_id_b)
                        scores = getInteractType(link.evidence_scores)
                        G.add_edge(i.node_id_b, a.node_id_b, score=link.combined_score,
                                   Homology=scores["Homology"], Experiment=scores["Experiment"],
                                   Database=scores["Database"], Textmining=scores["Textmining"],
                                   Genfusion=scores["Genfusion"], Coocurence=scores["Coocurence"],
                                   Neighborhood=scores["Neighborhood"], Coexpression=scores["Coexpression"], clicked=1, checkbox=1)
                        maxScore = link.combined_score if maxScore < link.combined_score else maxScore
                        minScore = link.combined_score if minScore > link.combined_score else minScore
                    except ObjectDoesNotExist:
                        pass

    d = json_graph.node_link_data(G)
    json_file = json.dumps(d)
    json_file = calcColorRange(json_file, maxhood, minhood, maxScore, minScore)
    return json_file

def calcColorRange(jsonFile, maxHood, minHood, maxScore, minScore):
    colorrangeNode = (maxHood-minHood)/4
    colorrangeLink = (maxScore-minScore)/4

    linkcolor = '"linkcolor": [{"color": "red", "value": '+str(minScore)+'}, ' \
                                '{"color": "yellow", "value": '+str(minScore+colorrangeLink)+'}, ' \
                                '{"color": "green", "value": '+str(maxScore-colorrangeLink)+'}, ' \
                                '{"color": "navy", "value": '+str(maxScore)+'}]'
    nodecolor = '"nodecolor": [{"color": "black", "value": 0}, ' \
                                '{"color": "blue", "value": '+str(minHood)+'}, ' \
                                '{"color": "green", "value": '+str(minHood+colorrangeNode)+'}, ' \
                                '{"color": "yellow", "value": '+str(maxHood-colorrangeNode)+'}, ' \
                                '{"color": "red", "value": '+str(maxHood)+'}]'
    searchString = '"Experiment":'
    jsonFile = jsonFile.replace(searchString, '"filterscore": 0, '+searchString)
    jsonFile = jsonFile[:-1]+', '+linkcolor+', '+nodecolor+'}'
    return jsonFile




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
            G.add_node(prot.protein_id, name=prot.preferred_name, hood=0)
            print prot.preferred_name
            for i in interactions:
                linked = i.node_id_b
                G.add_node(linked, name=Proteins.objects.get(protein_id=linked).preferred_name, hood=len(
                    NodeNodeLinks.objects.filter(node_id_a=linked).order_by('-combined_score')))
                G.add_edge(prot.protein_id, linked, score=i.combined_score, opacity=1)

            ''' Looking at each Protein and connecting them, if they interact'''
            for i in interactions:
                for a in interactions:
                    if not a.node_id_b == i.node_id_b:
                        if not (G.has_edge(a.node_id_b, i.node_id_b)):
                            if NodeNodeLinks.objects.filter(node_id_a=i.node_id_b, node_id_b=a.node_id_b).count() != 0:
                                link = NodeNodeLinks.objects.get(node_id_a=i.node_id_b, node_id_b=a.node_id_b)
                                G.add_edge(i.node_id_b, a.node_id_b, score=link.combined_score, opacity=0.2)
            print "Schreibe Datei"
            d = json_graph.node_link_data(G)
            data = {"url": 'cyanofactory/stringdb/protein_network/' + str(prot.protein_id) + '.json'}
            json.dump(d, open(data["url"], 'w'))
            a = prot.preferred_name
        else:
            print prot.preferred_name + " already exists"
    a = len(proteins)
    return HttpResponse("Finished - " + str(a))
    #return chelpers.render_queryset_to_response(request, data=data, template="stringdb\layout.html")


'''
Comparing 2 Proteins and checking if they have equal interactions partners
Until now 2 Proteins are related if the half of there interactions partners are equal
To Do: Working on a weight function to include Score and amount of interaction partner per Protein
'''


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


def getwholenetwork(request):
    print "Running"
    proteins = Proteins.objects.filter(species_id=1148)
    G = nx.Graph()
    amount = len(proteins)
    current = 1
    i = 0
    for prot in proteins:
        if current == amount/4:
            i += current
            print str(i)+" done"
            current = 1
        #print prot.preferred_name
        protid = prot.protein_id
        if not G.has_node(protid):
            G.add_node(protid, name=prot.preferred_name)
            interacts = NodeNodeLinks.objects.filter(node_id_a=protid)
            for interact in interacts:
                interactprot = Proteins.objects.get(protein_id=interact.node_id_b)
                if not G.has_node(interactprot.protein_id):
                    G.add_node(interactprot.protein_id, name=interactprot.preferred_name)
                if not G.has_edge(interactprot, protid):
                    G.add_edge(protid, interactprot.protein_id, score=interact.combined_score)
        current += 1
    d = json_graph.node_link_data(G)
    json.dump(d, open('cyanofactory/stringdb/protein_network_v1/cyano.json', 'w'))
    print "finish"
    return HttpResponse("Done")


def listInteracts(request, proteinid):
    linked_proteins = NodeNodeLinks.objects.filter(node_a_id=proteinid).order_by("-combined_score")


def listInteractions(protID):
    interactions = NodeNodeLinks.objects.filter(node_id_a=protID).order_by("-combined_score")
    for interact in interactions:
        interactType = getInteractType(interact.evidence_scores)
        interact.preferred_name = Proteins.objects.get(protein_id=interact.node_id_b).preferred_name
        prot = Proteins.objects.get(protein_id=interact.node_id_b).annotation

        if prot.__contains__(";"):
            name, rest = prot.split(";", 1)
        else:
            name = prot
        interact.annotation = name
        interact.evidence = interactType
    return interactions


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
        #elif scoretype == 6 or scoretype == 7:
        elif scoretype == 6:
            interactType["CoexpressionA"] = score[1]
        elif scoretype == 7:
            interactType["CoexpressionB"] = score[1]
        #elif scoretype == 8 or scoretype == 9:
        elif scoretype == 8:
            interactType["ExperimentA"] = score[1]
        elif scoretype == 9:
            interactType["ExperimentB"] = score[1]
        #elif scoretype == 10 or scoretype == 11:
        elif scoretype == 10:
            interactType["DatabaseA"] = score[1]
        elif scoretype == 11:
            interactType["DatabaseB"] = score[1]
        #elif scoretype == 12 or scoretype == 13:
        elif scoretype == 12:
            interactType["TextminingA"] = score[1]
        elif scoretype == 13:
            interactType["TextminingB"] = score[1]
        #elif scoretype == 14 or scoretype == 1 or scoretype == 2:
        elif scoretype == 1:
            interactType["NeighborhoodA"] = score[1]
        elif scoretype == 2:
            interactType["NeighborhoodB"] = score[1]
        elif scoretype == 14:
            interactType["NeighborhoodC"] = score[1]
        #elif scoretype == 15 or scoretype == 3:
        elif scoretype == 3:
            interactType["GenfusionA"] = score[1]
        elif scoretype == 15:
            interactType["GenfusionB"] = score[1]
        #elif scoretype == 16 or scoretype == 4:
        elif scoretype == 4:
            interactType["CoocurenceA"] = score[1]
        elif scoretype == 16:
            interactType["CoocurenceB"] = score[1]

    interactType["Coexpression"] = (1-((1-interactType["CoexpressionA"]/1000.0)*(1-interactType["CoexpressionB"]/1000.0)))*1000.0
    interactType["Experiment"] = (1-((1-interactType["ExperimentA"]/1000.0)*(1-interactType["ExperimentB"]/1000.0)))*1000.0
    interactType["Database"] = (1-((1-interactType["DatabaseA"]/1000.0)*(1-interactType["DatabaseB"]/1000.0)))*1000.0
    interactType["Textmining"] = (1-((1-interactType["TextminingA"]/1000.0)*(1-interactType["TextminingB"]/1000.0)))*1000.0
    interactType["Neighborhood"] = (1-((1-interactType["NeighborhoodA"]/1000.0)*(1-interactType["NeighborhoodB"]/1000.0)*(1-interactType["NeighborhoodC"]/1000.0)))*1000.0
    interactType["Genfusion"] = (1-((1-interactType["GenfusionA"]/1000.0)*(1-interactType["GenfusionB"]/1000.0)))*1000.0
    interactType["Coocurence"] = (1-((1-interactType["CoocurenceA"]/1000.0)*(1-interactType["CoocurenceB"]/1000.0)))*1000.0

    return interactType

def determineDistance(jsonfile, protA, protB):
    jsondata = json.load(jsonfile)
    g = json_graph.node_link_graph(jsondata)
    nProtA = Proteins.objects.get(preferred_name=protA, species_id=1148)
    nProtB = Proteins.objects.get(preferred_name=protB, species_id=1148)
    result = nx.shortest_path(g, nProtA.protein_id, nProtB.id)
    print result