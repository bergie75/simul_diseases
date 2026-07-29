import numpy as np
from networkx.generators.random_graphs import erdos_renyi_graph, barabasi_albert_graph
import networkx as nx
import matplotlib.pyplot as plt
import scipy

# scale = 10**3
# num_nodes = 40*scale
# G=barabasi_albert_graph(num_nodes, int(0.3/scale*num_nodes))
# adjacency_matrix = nx.to_scipy_sparse_array(G)
# degree_values = np.sum(adjacency_matrix, axis=0)

# plt.hist(degree_values, bins=100)
# plt.show()

rng = np.random.default_rng()

def construct_adj_mat(num_nodes, num_clusters, inter_prob, ba_m):
    if num_nodes % num_clusters != 0:
        raise ValueError("Cannot divide clusters evenly")
    nodes_per_cluster = num_nodes //num_clusters

    # adjacency matrix to return
    adj_mat = np.zeros((num_nodes, num_nodes))

    for i in range(0, num_clusters):
        for j in range(0, num_clusters):
            if i==j:
                G = barabasi_albert_graph(nodes_per_cluster, ba_m)
                sub_matrix = nx.to_numpy_array(G)
                adj_mat[i*nodes_per_cluster:(i+1)*nodes_per_cluster, i*nodes_per_cluster:(i+1)*nodes_per_cluster] = sub_matrix
            else:
                sub_matrix = np.random.rand(nodes_per_cluster, nodes_per_cluster) <= inter_prob
                sub_trans = np.transpose(sub_matrix)
                adj_mat[i*nodes_per_cluster:(i+1)*nodes_per_cluster, j*nodes_per_cluster:(j+1)*nodes_per_cluster] = sub_matrix
                adj_mat[j*nodes_per_cluster:(j+1)*nodes_per_cluster, i*nodes_per_cluster:(i+1)*nodes_per_cluster] = sub_trans

    return adj_mat

if __name__ == "__main__":
    A = construct_adj_mat(40*10**3, 100, 0.0005, 30)
    degree_values = np.sum(A, axis=0)
    print(np.linalg.norm(A-np.transpose(A)))

    # plt.hist(degree_values, bins=100)
    # plt.show()