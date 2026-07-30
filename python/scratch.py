import numpy as np
import os
import networkx as nx
from networkx.generators.random_graphs import erdos_renyi_graph, newman_watts_strogatz_graph, barabasi_albert_graph, fast_gnp_random_graph

a=np.array([1,2,3,4])
b=np.array([5,4,3,2])

print([x>=b[j] for j,x in enumerate(a)])

        #G = newman_watts_strogatz_graph(num_nodes, int(num_nodes/100), 0)
        #G = erdos_renyi_graph(num_nodes, 0.7/scale)
        # G=barabasi_albert_graph(num_nodes, int(0.3/scale*num_nodes))
        # adjacency_matrix = nx.to_scipy_sparse_array(G)

print(os.path.join(os.getcwd(), "", "hello"))