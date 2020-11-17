import networkx as nx
import matplotlib.pyplot as plt
import random

n = 100
pos = {i: (random.uniform(0, 1), random.uniform(0, 1)) for i in range(n)}

node_list = [i for i in range(n)]

H = nx.Graph()
H.add_nodes_from(node_list, pos=pos)

# nx.draw(H, pos=pos, node_size=30)
# plt.show()

G = nx.random_geometric_graph(n, 0.2, pos=pos)

deg = [G.degree[i] for i in range(n)]
vmin = min(deg)
vmax = max(deg)
print(deg)

nx.draw(G, pos=pos, node_color=deg, node_size=30, width=0.5, cmap=plt.cm.Reds)
sm = plt.cm.ScalarMappable(cmap=plt.cm.Reds, norm=plt.Normalize(vmin = vmin, vmax=vmax))
plt.colorbar(sm)
plt.show()
