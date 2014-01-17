from django.http import HttpResponse
from django.shortcuts import render_to_response
from stringdb.models import Proteins, NodeNodeLinks
import networkx as nx
import json
from networkx.readwrite import json_graph
import cyano.helpers as chelpers
import html_server
# Create your views here.

def index(request):
    protein_list = Proteins.objects.order_by('-protein_id')[:50]
    proteincount = Proteins.objects.count()
    count = proteincount
    return HttpResponse(count)


def graph(request):
    proteins = Proteins.objects.filter(species_id=1148)
    print len(proteins)
    for prot in proteins:
        print "Lade Protein"
        print "Erstelle Graph"
        G = nx.Graph()
        interactions = NodeNodeLinks.objects.filter(node_id_a=prot.protein_id).order_by('-combined_score')[:20]
        #interactions = NodeNodeLinks.objects.filter(node_id_a=prot.protein_id)
        G.add_node(prot.protein_id, name=prot.preferred_name, score=1000)
        print prot.preferred_name
        print "Suche nach Interaktionspartnern"
        for i in interactions:
            print "ID - {}".format(i.node_id_b)
            linked = i.node_id_b
            G.add_node(linked, name=Proteins.objects.get(protein_id=linked).preferred_name, score=i.combined_score,
                       shape="circle")
            G.add_edge(prot.protein_id, linked, score=i.combined_score)
        for i in interactions:
            print "ID 1 - {}".format(i.node_id_b)
            for a in interactions:
                #print "ID 2 - {}".format(a.node_id_b)
                if not a.node_id_b == i.node_id_b:
                    #if not (NodeNodeLinks.objects.get(node_id_a=i.node_id_b, node_id_b=a.node_id_b).count() == 0):
                    if not (G.has_edge(a.node_id_b, i.node_id_b)):
                        if NodeNodeLinks.objects.filter(node_id_a=i.node_id_b, node_id_b=a.node_id_b).count() != 0:
                            link = NodeNodeLinks.objects.get(node_id_a=i.node_id_b, node_id_b=a.node_id_b)
                            G.add_edge(i.node_id_b, a.node_id_b, score=link.combined_score)
        print "Schreibe Datei"
        d = json_graph.node_link_data(G)
        data = {"url": 'cyanofactory/stringdb/network/' + prot.preferred_name + '.json'}
        json.dump(d, open(data["url"], 'w'))
        a = prot.preferred_name
    return HttpResponse("Finished")
    #return chelpers.render_queryset_to_response(request, data=data, template="stringdb\layout.html")



def compareInteraction(request):
    proteins = Proteins.objects.filter(species_id=1148).protein_id
    connected_prots = {}
    for prot_a in proteins:
        connected = []
        counter = 0
        linked_a = NodeNodeLinks.objects.filter(node_a_id=prot_a).order_by('-combined_score').node_b_id[:10]
        for prot_b in proteins:
            if not (prot_a == prot_b):
                linked_b = NodeNodeLinks.objects.filter(node_a_id=prot_b).order_by('-combined_score').node_b_id[:10]
                for link in linked_a:
                    if (link in linked_b):
                        counter += 1
                if (linked_a / 2 <= counter):
                    connected.append(prot_b)
        if not (0 == len(connected)):
            connected_prots[prot_a] = connected
    return HttpResponse("Done")


def getwholenetwork(request):
    proteins = Proteins.objects.filter(species_id=1148)[:2]
    G = nx.Graph()
    for prot in proteins:
        protid = prot.protein_id
        if not (G.has_node(protid)):
            G.add_node(protid, name=prot.preferred_name)
            interacts = NodeNodeLinks.objects.filter(node_a_id=protid)
            for interact in interacts:
                interactprot = Proteins.objects.get(protein_id=interact.node_id_b)
                if not (G.has_node(interactprot.protein_id)):
                    G.add_node(interactprot.protein_id, name=interactprot.preferred_name)
                if not (G.has_edge(interactprot, protid)):
                    G.add_edge(protid, interactprot.protein_id, score=interact.combined_score)
    d = json_graph.node_link_data(G)
    json.dump(d, open('cyanofactory/stringdb/network/cyano.json', 'w'))


