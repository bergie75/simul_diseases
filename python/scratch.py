import numpy as np
import os
import networkx as nx
from networkx.generators.random_graphs import erdos_renyi_graph, newman_watts_strogatz_graph, barabasi_albert_graph, fast_gnp_random_graph

name_mod = "sine_two_week"
cwd = os.getcwd()
savefile = os.path.join(cwd, "parameters", f"{name_mod}.npy")

num_days = 90
trans_mods = np.array([0.5+0.5*np.cos(2*np.pi*t/14) for t in range(0, 90)])
np.save(savefile, trans_mods)