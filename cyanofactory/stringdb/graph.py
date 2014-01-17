import networkx as nx
import json
from networkx.readwrite import json_graph
import html_server
__author__ = 'Eric'

G = nx.barbell_graph(6, 3)
for n in G:
    G.node[n]['name'] = n
d = json_graph.node_link_data(G)
json.dump(d, open('force.json', 'w'))
html_server.load_url('force.json')
