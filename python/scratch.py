import numpy as np
import networkx as nx
from networkx.generators.random_graphs import erdos_renyi_graph, newman_watts_strogatz_graph, barabasi_albert_graph, fast_gnp_random_graph

a=np.array([1,2,3,4])
b=np.array([5,4,3,2])

print([x>=b[j] for j,x in enumerate(a)])