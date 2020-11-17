import networkx as nx
import matplotlib.pyplot as plt
import random

n = 30
p = 0.2

node_list = [i for i in range(n)]
edge_list = [(i,j) for i in range(n) for j in range(i+1, n) if random.random()<p]

G = nx.Graph()

G.add_nodes_from(node_list)
G.add_edges_from(edge_list)

nx.draw(G, node_color='darkred', width=0.3, node_size=70)

print(len(G.nodes()))

plt.show()