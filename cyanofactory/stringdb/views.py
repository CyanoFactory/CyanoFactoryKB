from django.http import HttpResponse
from django.shortcuts import render_to_response
from stringdb.models import Proteins, NodeNodeLinks, ProteinsOrthgroups
import networkx as nx
import json
from networkx.readwrite import json_graph
import os
from cyano.models import ProteinComparison
from cyano.helpers import render_queryset_to_response

# Create your views here.

def checkInteraction(request, id):
   # json_file = proteingraph(id)
    interactions = listInteractions(id)
    return render_queryset_to_response(
        request=request,
        template="stringdb/index.html",
        data={
            'protid': id,
            'prot_name': Proteins.objects.get(protein_id=id).preferred_name,
            'interacts': interactions
            #'json': json_file
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


def proteingraph(protein_id):
    print protein_id
    protein = Proteins.objects.get(protein_id=protein_id)
    G = nx.Graph()
    interactions = NodeNodeLinks.objects.filter(node_id_a=protein.protein_id).order_by('-combined_score')[:10]
    #interactions = NodeNodeLinks.objects.filter(node_id_a=prot.protein_id)
    G.add_node(protein.protein_id, name=protein.preferred_name, hood=0)
    print protein.preferred_name
    print "Anzahl"+str(len(interactions))
    print "Suche nach Interaktionspartnern"
    for i in interactions:
        # print "ID - {}".format(i.node_id_b)
        print "Betrachte Partner"
        linked = i.node_id_b
        G.add_node(linked, name=Proteins.objects.get(protein_id=linked).preferred_name, hood=len(
            NodeNodeLinks.objects.filter(node_id_a=linked).order_by('-combined_score')))
        G.add_edge(protein.protein_id, linked, score=i.combined_score, opacity=1)
    for i in interactions:
        # print "ID 1 - {}".format(i.node_id_b)
        for a in interactions:
            #print "ID 2 - {}".format(a.node_id_b)
            if not a.node_id_b == i.node_id_b:
                #if not (NodeNodeLinks.objects.get(node_id_a=i.node_id_b, node_id_b=a.node_id_b).count() == 0):
                if not (G.has_edge(a.node_id_b, i.node_id_b)):
                    if NodeNodeLinks.objects.filter(node_id_a=i.node_id_b, node_id_b=a.node_id_b).count() != 0:
                        link = NodeNodeLinks.objects.get(node_id_a=i.node_id_b, node_id_b=a.node_id_b)
                        G.add_edge(i.node_id_b, a.node_id_b, score=link.combined_score, opacity=0.2)
    print "Schreibe Datei"
    d = json_graph.node_link_data(G)
    json_file = json.dump(d)
    return json_file


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
            #print "Lade Protein"
            #print "Erstelle Graph"
            print str(counter) + " - " + str(len(proteins))
            G = nx.Graph()
            interactions = NodeNodeLinks.objects.filter(node_id_a=prot.protein_id).order_by('-combined_score')[:10]
            #interactions = NodeNodeLinks.objects.filter(node_id_a=prot.protein_id)
            G.add_node(prot.protein_id, name=prot.preferred_name, hood=0)
            print prot.preferred_name
            #print "Suche nach Interaktionspartnern"
            for i in interactions:
                # print "ID - {}".format(i.node_id_b)
                linked = i.node_id_b
                G.add_node(linked, name=Proteins.objects.get(protein_id=linked).preferred_name, hood=len(
                    NodeNodeLinks.objects.filter(node_id_a=linked).order_by('-combined_score')))
                G.add_edge(prot.protein_id, linked, score=i.combined_score, opacity=1)
            for i in interactions:
                # print "ID 1 - {}".format(i.node_id_b)
                for a in interactions:
                    #print "ID 2 - {}".format(a.node_id_b)
                    if not a.node_id_b == i.node_id_b:
                        #if not (NodeNodeLinks.objects.get(node_id_a=i.node_id_b, node_id_b=a.node_id_b).count() == 0):
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
    proteins = Proteins.objects.filter(species_id=1148)
    G = nx.Graph()
    for prot in proteins:
        print prot.preferred_name
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
    d = json_graph.node_link_data(G)
    json.dump(d, open('cyanofactory/stringdb/protein_network/cyano.json', 'w'))
    return HttpResponse("Done")


def listInteracts(request, proteinid):
    linked_proteins = NodeNodeLinks.objects.filter(node_a_id=proteinid).order_by("-combined_score")

def listInteractions(protID):
    interactions = NodeNodeLinks.objects.filter(node_id_a=protID).order_by("-combined_score")
    for interact in interactions:
        interact.preferred_name = Proteins.objects.get(protein_id=interact.node_id_b).preferred_name
    return interactions
