import numpy as np
import os
from networkx.generators.random_graphs import erdos_renyi_graph, newman_watts_strogatz_graph, barabasi_albert_graph
import matplotlib.pyplot as plt
import networkx as nx
from node import Node
from scipy.io import mmread

rng = np.random.default_rng()
status_to_num = {"S": 0, "E": 1, "I": 2, "R": 3}

def build_node_list(adjacency_matrix, starting_exposure_frac,
                    transmit_prob, fall_ill_prob, recover_prob,
                     dis_weight,prior_res):
    
    # check for square adjacency matrix
    m,n = adjacency_matrix.shape
    if m != n:
        raise ValueError("Dimensional mismatch, invalid matrix")
    
    node_list = []
    all_exposures = [np.zeros(n),np.zeros(n)]
    for i in range(0,n):
        # randomly choose a disease to update first so that we are not biased
        first_disease = int((rng.uniform() <= 0.5))  # is either 0 or 1 (as Boolean) with 50% probability
        
        generated_status = ["S", "S"]  # default is exposed to neither disease, only get exposed to one initially
        blocking_disease = None
        if rng.uniform() <= starting_exposure_frac[first_disease]:
            generated_status[first_disease] = "E"
            blocking_disease = first_disease
            all_exposures[first_disease][i] = 1  # keep track of who is initially infected
        elif rng.uniform() <= starting_exposure_frac[1-first_disease]:
            generated_status[1-first_disease] = "E"
            blocking_disease = 1-first_disease
            all_exposures[1-first_disease][i] = 1  # keep track of who is initially infected
        
        new_node = Node(i, transmit_prob=transmit_prob, fall_ill_prob=fall_ill_prob, 
                        recover_prob=recover_prob, status=generated_status, blocking_disease=blocking_disease,
                        inf_decision_weights=dis_weight, prior_resistance=prior_res)
        node_list.append(new_node)
    
    # after all nodes are generated, we give them their neighbors from the adjacency list
    for Current_node in node_list:
        neighbor_list = []
        for n_index in adjacency_matrix[Current_node.index,:].nonzero()[0]:
            neighbor_list.append(node_list[n_index])
        Current_node.update_neighbors(neighbor_list)
    
    return node_list, all_exposures

def late_exposure(node_list, disease_index, exposure_frac):
    for node in node_list:
        if node.blocking_disease is None and rng.uniform() <= exposure_frac:
            node.blocking_disease = disease_index
            node.status[disease_index] = "E"

if __name__ == "__main__":
    # parameter folder
    cwd = os.getcwd()
    par_folder = os.path.join(cwd, "parameters")
    new_matrix = True

    # will keep track of all daily statistics
    daily_counts = [[], []]
    affected_nodes = [set(), set()]
    global_blocking_events = []
    count_self_blocks = True
    
    # other parameters
    num_days = 150
    scale = 10**3
    num_nodes = 40*scale
    total_node_set = set(range(0,num_nodes))
    starting_exposure_frac = np.array([0.0, 1.0])/scale
    delay_day = 0

    # delayed entry of second disease options
    if delay_day > 0:
        starting_exposure_frac[-1] = 0
    late_frac = 1.0/scale

    num_diseases = len(starting_exposure_frac)
    trans_prob = [0.04, 0.01]
    fall_ill = [0.4,0.1]  # prob of becoming fully infected
    rec_prob = [2.0/num_days,2.0/num_days]
    dis_weight = [1,1]
    prior_res = [0,0.4]

    # adjacency matrix construction and/or saving
    if new_matrix:
        #G = newman_watts_strogatz_graph(num_nodes, int(num_nodes/100), 0)
        #G = erdos_renyi_graph(num_nodes, 0.7/scale)
        G=barabasi_albert_graph(num_nodes, int(0.7/scale*num_nodes))
        adjacency_matrix = nx.to_scipy_sparse_array(G)
        np.save(os.path.join(par_folder, "adj_mat"), adjacency_matrix)
        
    else:
        adjacency_matrix = np.load(os.path.join(par_folder, "adj_mat.npy"))

    node_list, all_exposures = build_node_list(adjacency_matrix, starting_exposure_frac,
                                                trans_prob, fall_ill, rec_prob, dis_weight, prior_res)

    disease_one_exposed = np.sum(all_exposures[0])
    disease_one_initial = np.array([num_nodes-disease_one_exposed,disease_one_exposed,0,0])
    disease_two_exposed = np.sum(all_exposures[1])
    disease_two_initial = np.array([num_nodes-disease_two_exposed,disease_two_exposed,0,0])

    daily_counts[0].append(disease_one_initial)
    daily_counts[1].append(disease_two_initial)

    for day in range(0, num_days):
        print(f"Day {1+day}")
        if delay_day > 0 and delay_day-1 == day:
            late_exposure(node_list, 1, late_frac)
            print("Late exposure to second disease initiated.")
        
        # progress the status of the various nodes, and transmit disease among neighbors
        for Current_node in node_list:
            Current_node.progress_status()
        
        # transmission events
        for Current_node in node_list:
            for disease_index in range(0, num_diseases):
                blocked_trans = Current_node.transmit(disease_index, day=day)
                global_blocking_events.extend(blocked_trans)
        
        # end of day update of status
        current_status_totals = [np.array([0,0,0,0]), np.array([0,0,0,0])]
        
        # choose which transmission occured if multiple, and collect stats
        for Current_node in node_list:
            chosen_disease = Current_node.choose_exposure_event()
            if chosen_disease is not None:
                affected_nodes[chosen_disease].add(Current_node.index)
            for disease_index in range(0, num_diseases):
                current_status_totals[disease_index][status_to_num[Current_node.status[disease_index]]] += 1
        
        # add aggregated status totals to list to track disease over time
        daily_counts[0].append(current_status_totals[0])
        daily_counts[1].append(current_status_totals[1])

    # extract case counts in plottable format
    disease_one_counts = np.array(daily_counts[0])
    disease_two_counts = np.array(daily_counts[1])

    # draw graph
    plot_graph = False
    if plot_graph:
        pos = nx.spring_layout(G, seed=3113794652)
        nx.draw_networkx_nodes(G, pos, nodelist=list(affected_nodes[0]), node_color="blue")
        nx.draw_networkx_nodes(G, pos, nodelist=list(affected_nodes[1]), node_color="orange")
        remaining_nodes = total_node_set.difference(affected_nodes[0].union(affected_nodes[1]))
        nx.draw_networkx_nodes(G, pos, nodelist=list(remaining_nodes), node_color="gray")
        nx.draw_networkx_edges(G, pos, width=0.3, alpha=0.5)
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)
    
    ax1.plot(disease_one_counts[:,1], label="Disease one exposed")
    ax1.plot(disease_two_counts[:,1], label="Disease two exposed")
    ax2.plot(disease_one_counts[:,2], label="Disease one infected")
    ax2.plot(disease_two_counts[:,2], label="Disease two infected")
    ax3.plot(disease_one_counts[:,3], label="Disease one recovered")
    ax3.plot(disease_two_counts[:,3], label="Disease two recovered")
    ax4.plot(disease_one_counts[:,0], label="Disease one susceptible")
    ax4.plot(disease_two_counts[:,0], label="Disease two susceptible")    
    ax1.legend()
    ax2.legend()
    ax3.legend()
    ax4.legend()
    plt.show()

    fig2, ((pi1, pi2), (hist, bar)) = plt.subplots(2,2)
    exposure_blocks = np.array([0,0])
    infection_blocks = np.array([0,0])
    block_day_histogram = np.zeros((2,num_days))
    nodes_issuing_blocks = np.zeros(num_nodes)

    # post-processing on blocking events
    for status_readout in global_blocking_events:
        block_issued, disease_index, index, status, day, self_block = status_readout
        if block_issued and (count_self_blocks or not self_block):
            block_day_histogram[disease_index, day] += 1
            if status[disease_index] == "E":
                exposure_blocks[disease_index] += 1
            elif status[disease_index] == "I":
                infection_blocks[disease_index] += 1
    
    disease_one_blocks = exposure_blocks[0] + infection_blocks[0]
    if disease_one_blocks > 0:
        disease_one_block_fracs = [exposure_blocks[0]/disease_one_blocks, infection_blocks[0]/disease_one_blocks]
        pi1.pie(disease_one_block_fracs, labels=["Exposures", "Infections"], colors= ["red", "black"])
        pi1.legend(loc='center left', bbox_to_anchor=(1.05, 0.9))
    
    disease_two_blocks = exposure_blocks[1] + infection_blocks[1]
    if disease_two_blocks > 0:
        disease_two_block_fracs = [exposure_blocks[1]/disease_two_blocks, infection_blocks[1]/disease_one_blocks]
        pi2.pie(disease_two_block_fracs, labels = ["Exposures", "Infections"], colors= ["red", "black"])
        pi2.legend(loc='center left', bbox_to_anchor=(1.05, 0.9))

    hist.stairs(block_day_histogram[0,:], fill=True, label="Disease 1", alpha=0.3)
    hist.stairs(block_day_histogram[1,:], fill=True, label="Disease 2", alpha=0.3)
    
    bar.bar([-1,1], [disease_one_blocks, disease_two_blocks], tick_label=["Disease 1", "Disease 2"], color= ["blue", "orange"])

    plt.show()