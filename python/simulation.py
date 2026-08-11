import numpy as np
import os
from sys import exit as leave_program
import yaml
import pickle
from copy import deepcopy
from networkx.generators.random_graphs import barabasi_albert_graph
import matplotlib.pyplot as plt
import networkx as nx
from node import Node
from matplotlib.transforms import ScaledTranslation

rng = np.random.default_rng()
status_to_num = {"S": 0, "I": 1, "H": 2, "R": 3}

def modded_division(x,y):
    result = []
    for i in range(0, len(x)):
        if x[i] == 0:
            result.append(0)
        else:
            result.append(x[i]/y[i])

    return result

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

def neighbor_only_node_list(adjacency_matrix):
    # check for square adjacency matrix
    m,n = adjacency_matrix.shape
    if m != n:
        raise ValueError("Dimensional mismatch, invalid matrix")
    
    node_list = []
    for i in range(0,n):
        new_node = Node(i)  # should use generic values for all parameters from the config list
        node_list.append(new_node)
    
    # after all nodes are generated, we give them their neighbors from the adjacency list
    for Current_node in node_list:
        neighbor_list = []
        for n_index in adjacency_matrix[Current_node.index,:].nonzero()[0]:
            neighbor_list.append(n_index)
        Current_node.update_neighbors(neighbor_list)
    
    return node_list

def reset_nodes_new_sim(node_list, starting_exposure_frac, config_opts={}):

    # track information on initial infections
    n=len(node_list)
    all_exposures = [np.zeros(n),np.zeros(n)]

    for i,Current_node in enumerate(node_list):
        # modify any variables which must be altered from default values
        Current_node.dict_update(config_opts)

        # wipe any effects from previous simulations
        Current_node.set_status(0, "S")
        Current_node.set_status(1, "S")
        Current_node.set_block(None)

        # randomly choose a disease to update first so that we are not biased
        first_disease = int((rng.uniform() <= 0.5))  # is either 0 or 1 (as Boolean) with 50% probability

        if rng.uniform() <= starting_exposure_frac[first_disease]:
            Current_node.set_status(first_disease, "I")
            Current_node.set_block(first_disease)
            all_exposures[first_disease][i] = 1  # keep track of who is initially infected
        elif rng.uniform() <= starting_exposure_frac[1-first_disease]:
            Current_node.set_status(1-first_disease, "I")
            Current_node.set_block(1-first_disease)
            all_exposures[1-first_disease][i] = 1  # keep track of who is initially infected
    
    return all_exposures

def late_exposure(node_list, disease_index, exposure_frac):
    for node in node_list:
        if node.blocking_disease is None and rng.uniform() <= exposure_frac:
            node.blocking_disease = disease_index
            node.status[disease_index] = "I"

def run_simulation(run_name, opt_args={}, node_folder="", trans_mod=None):
    cwd = os.getcwd()
    save_folder = os.path.join(cwd, "output", run_name)
    par_folder = os.path.join(cwd, "parameters")

    print(f"Working on {run_name}")

    # locate configuration file
    with open(os.path.join(par_folder, "config.yaml"), "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    # after loading config file, use mandatory changes
    for opt_arg in opt_args.keys():
        cfg[opt_arg] = deepcopy(opt_args[opt_arg])

    num_days = cfg["num_days"]
    scale = cfg["scale"]
    delay_day = cfg["delay_day"]
    new_matrix = cfg["new_matrix"]
    count_self_blocks = cfg["count_self_blocks"]

    # allow for weather and other modifications of transmission probability
    if trans_mod is not None:
        trans_mod_file = os.path.join(cwd, "parameters", f"{trans_mod}".npy)
        trans_modifiers = np.load(trans_mod_file)
        cfg["trans_mod"] = trans_mod

        # check for compatibility with num_days
        if len(trans_modifiers) != num_days:
            raise ValueError("Transmission modifiers does not occur over the appropriate number of days")
    else:
        trans_modifiers = np.ones(num_days)

    # for convenience so node_folder does not need to be specified when creating a new matrix
    if new_matrix and node_folder == "":
        node_folder = run_name

    # compute variables after loading config values
    num_nodes = 40*scale
    starting_exposure_frac = np.array([1.0, 1.0])/scale

    # will keep track of all daily statistics
    daily_counts = [[], []]
    affected_nodes = [set(), set()]
    global_blocking_events = []
    used_recovered_node = np.zeros((2,num_days))
    used_recovered_degree = np.zeros((2,num_days))

    # delayed entry of second disease options
    if delay_day > 0:
        starting_exposure_frac[-1] = 0
    late_frac = 1.0/scale

    # define transmission parameters
    num_diseases = len(starting_exposure_frac)
    rec_weights = cfg["rec_weights"]
    rec_prob = [x/num_days for x in rec_weights]

    # this is to get around a weird issue with recovery probability
    sim_cfg = {}
    node_attributes = ["dis_weight", "fall_ill", "prior_res", "trans_prob", "large_resp"]  # only keep these variables to pass to sim_cfg
    for attribute in node_attributes:
        sim_cfg[attribute] = deepcopy(cfg[attribute])
    sim_cfg["rec_prob"] = rec_prob  # this is what is needed by the nodes

    # matrix parameters
    m_BA = cfg["m_BA"]
    scaled_p_ER = cfg["scaled_p_ER"]
    n_clusters = cfg["n_clusters"]

    # adjacency matrix construction and/or saving
    print("Generating network information ...")
    if new_matrix:
        if not os.path.isdir(os.path.join(par_folder, node_folder)):
            os.mkdir(os.path.join(par_folder, node_folder))
        else:
            response = input("Directory already exists. Overwrite? (y/n): ")
            if response != "y":
                leave_program()

        adjacency_matrix = construct_adj_mat(num_nodes, n_clusters, scaled_p_ER/scale, m_BA)
        degree_values = np.sum(adjacency_matrix, axis=0)
        np.save(os.path.join(par_folder, node_folder, "degree_values"), degree_values)

        node_list = neighbor_only_node_list(adjacency_matrix)
        with open(os.path.join(par_folder, node_folder, "node_list"), 'wb') as f:
            pickle.dump(node_list, f)

    else:
        degree_values = np.load(os.path.join(par_folder, node_folder, "degree_values.npy"))
        with open(os.path.join(par_folder, node_folder, "node_list"), 'rb') as f:
            node_list = pickle.load(f)

    # output degree information
    print(f"Maximum degree: {max(degree_values)}")
    print(f"Mean degree value: {np.mean(degree_values)}\n")

    print("Generating initial infections ...\n")
    all_exposures = reset_nodes_new_sim(node_list, starting_exposure_frac, config_opts=sim_cfg)  # enforces all parameter value changes

    disease_one_exposed = np.sum(all_exposures[0])
    disease_one_initial = np.array([num_nodes-disease_one_exposed,disease_one_exposed,0,0])
    disease_two_exposed = np.sum(all_exposures[1])
    disease_two_initial = np.array([num_nodes-disease_two_exposed,disease_two_exposed,0,0])

    daily_counts[0].append(disease_one_initial)
    daily_counts[1].append(disease_two_initial)

    for day in range(0, num_days):
        print(f"Day {1+day}")
        # determine impact of exogeneous factors on transmission
        daily_trans_mod = trans_modifiers[day]

        if delay_day > 0 and delay_day-1 == day:
            late_exposure(node_list, 1, late_frac)
            print("Late exposure to second disease initiated.")
        
        # progress the status of the various nodes, and transmit disease among neighbors
        for Current_node in node_list:
            Current_node.progress_status()
        
        # transmission events
        for Current_node in node_list:
            for disease_index in range(0, num_diseases):
                blocked_trans = Current_node.transmit(node_list, disease_index, day=day, trans_mod=daily_trans_mod)
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

    # post-processing after the run on blocking events
    exposure_blocks = np.array([0,0])
    infection_blocks = np.array([0,0])
    block_day_histogram = np.zeros((2,num_days))

    for status_readout in global_blocking_events:
        block_issued, disease_index, index, status, day, self_block = status_readout
        if block_issued and (count_self_blocks or not self_block):
            block_day_histogram[disease_index, day] += 1
            if status[1-disease_index] == "I":
                exposure_blocks[disease_index] += 1
            elif status[1-disease_index] == "H":
                infection_blocks[disease_index] += 1

        if not block_issued and status[1-disease_index] == "R":
            used_recovered_node[disease_index, day] += 1
            used_recovered_degree[disease_index, day] += degree_values[index]

    # save all results
    if not os.path.exists(save_folder):
        os.mkdir(save_folder)

    # numpy save in newly-created directory
    np.save(os.path.join(save_folder, "disease_one_counts.npy"), disease_one_counts)
    np.save(os.path.join(save_folder, "disease_two_counts.npy"), disease_two_counts)
    np.save(os.path.join(save_folder, "exposure_blocks.npy"), exposure_blocks)
    np.save(os.path.join(save_folder, "infection_blocks.npy"), infection_blocks)
    np.save(os.path.join(save_folder, "used_recovered_node.npy"), used_recovered_node)
    np.save(os.path.join(save_folder, "used_recovered_degree.npy"), used_recovered_degree)
    np.save(os.path.join(save_folder, "block_day_histogram.npy"), block_day_histogram)
    np.save(os.path.join(save_folder, "degree_values.npy"), degree_values)

    # add a note about which node folder was used, since graph parameters may not be changed when using a saved network
    cfg["node_folder"] = node_folder

    # save copy of config file
    with open(os.path.join(save_folder, "config_copy.yaml"), 'w') as f:
        yaml.dump(cfg, f)

def re_run_simulation(run_name):
    cwd = os.getcwd()
    par_folder = os.path.join(cwd, "output", run_name)

    # locate configuration file
    with open(os.path.join(par_folder, "config_copy.yaml"), "r", encoding="utf-8") as f:
        opt_args = yaml.safe_load(f)

    node_folder = opt_args["node_folder"]
    opt_args["new_matrix"] = False

    # check for any modifications to transmission
    if "trans_mod" in opt_args.keys():
        trans_mod = opt_args["trans_mod"]
    else:
        trans_mod = None

    run_simulation(run_name, opt_args=opt_args, node_folder=node_folder, trans_mod=trans_mod)

def plot_outputs(run_name):
    cwd = os.getcwd()
    save_folder = os.path.join(cwd, "output", run_name)

    disease_one_counts = np.load(os.path.join(save_folder, "disease_one_counts.npy"))
    disease_two_counts = np.load(os.path.join(save_folder, "disease_two_counts.npy"))
    exposure_blocks = np.load(os.path.join(save_folder, "exposure_blocks.npy"))
    infection_blocks = np.load(os.path.join(save_folder, "infection_blocks.npy"))
    used_recovered_node = np.load(os.path.join(save_folder, "used_recovered_node.npy"))
    used_recovered_degree = np.load(os.path.join(save_folder, "used_recovered_degree.npy"))
    block_day_histogram = np.load(os.path.join(save_folder, "block_day_histogram.npy"))

    # create subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)
    fig2, ((pi1, pi2), (hist, bar)) = plt.subplots(2,2)

    ax1.plot(disease_one_counts[:,1], label="Disease 0 infected")
    ax1.plot(disease_two_counts[:,1], label="Disease 1 infected")
    # ax2.plot(disease_one_counts[:,2], label="Disease 0 infected")
    # ax2.plot(disease_two_counts[:,2], label="Disease 1 infected")
    ax2.plot(disease_one_counts[:,3], label="Disease 0 recovered")
    ax2.plot(disease_two_counts[:,3], label="Disease 1 recovered")
    ax3.plot(used_recovered_node[0,:], label="Disease 0 recovered nodes used")
    ax3.plot(used_recovered_node[1,:], label="Disease 1 recovered nodes used")
    # ax4.plot(used_recovered_degree[0,:], label="Disease 0 recovered degree used")
    # ax4.plot(used_recovered_degree[1,:], label="Disease 1 recovered degree used")
    ax4.plot(modded_division(used_recovered_degree[0,:], used_recovered_node[0,:]), label="Disease 0 recovered degree used")
    ax4.plot(modded_division(used_recovered_degree[1,:], used_recovered_node[1,:]), label="Disease 1 recovered degree used")
    for ax in (ax1, ax2, ax3, ax4):
        ax.grid(True)
        ax.legend()

    disease_one_blocks = exposure_blocks[0] + infection_blocks[0]
    if disease_one_blocks > 0:
        disease_one_block_fracs = [exposure_blocks[0]/disease_one_blocks, infection_blocks[0]/disease_one_blocks]
        pi1.pie(disease_one_block_fracs, labels=["Infections", "Heightened immunity"], colors= ["red", "black"])
        pi1.legend(loc='center left', bbox_to_anchor=(1.05, 0.9))
    
    disease_two_blocks = exposure_blocks[1] + infection_blocks[1]
    if disease_two_blocks > 0:
        disease_two_block_fracs = [exposure_blocks[1]/disease_two_blocks, infection_blocks[1]/disease_one_blocks]
        pi2.pie(disease_two_block_fracs, labels = ["Infections", "Heightened immunity"], colors= ["red", "black"])
        pi2.legend(loc='center left', bbox_to_anchor=(1.05, 0.9))

    hist.stairs(block_day_histogram[0,:], fill=True, label="Disease 0", alpha=0.3)
    hist.stairs(block_day_histogram[1,:], fill=True, label="Disease 1", alpha=0.3)
    hist.set_ylabel("Number of times disease was blocked")
    
    bar.bar([-1,1], [disease_one_blocks, disease_two_blocks], tick_label=["Disease 0", "Disease 1"], color= ["blue", "orange"])
    bar.set_ylabel("Number of times disease was blocked")

    plt.show()

def comparative_outputs(run_names, standardize_y_axis=False):
    # administrative variables
    num_to_plot = len(run_names)
    cwd = os.getcwd()
    caption = "A comparative analysis of "

    if num_to_plot == 4:
        letters = ["a)", "b)", "c)", "d)"]
        fig, ax = plt.subplots(2,2)
        peak_value = -np.inf

        for i, run_name in enumerate(run_names):
            row_num = i // 2
            col_num = i % 2

            # load results
            save_folder = os.path.join(cwd, "output", run_name)
            disease_one_counts = np.load(os.path.join(save_folder, "disease_one_counts.npy"))
            disease_two_counts = np.load(os.path.join(save_folder, "disease_two_counts.npy"))
            ax[row_num, col_num].plot(disease_one_counts[:,1], label="Disease 1")
            ax[row_num, col_num].plot(disease_two_counts[:,1], label="Disease 2")
            ax[row_num, col_num].grid(True)
            ax[row_num, col_num].legend()

            peak_value = max(peak_value, np.max(disease_one_counts[:,1]), np.max(disease_two_counts[:,2]))

            if col_num == 0:
                ax[row_num, col_num].set_ylabel("Infected nodes")

            if row_num == 1:
                ax[row_num, col_num].set_xlabel("Day")

            # add letterin
            ax[row_num, col_num].text(
                    -0.1, 1.15, letters[i], transform=(
                        ax[row_num, col_num].transAxes + ScaledTranslation(-25/72, -25/72, fig.dpi_scale_trans)),
                    fontsize='large', va='bottom', fontfamily='serif')

            # adjust caption
            if i<3:
                caption += f"{letters[i]} {run_names[i]}, "
            else:
                caption += f"and {letters[i]} {run_names[i]}."

        if standardize_y_axis:
                    for i in range(0, 4):
                        row_num = i // 2
                        col_num = i % 2
                        ax[row_num, col_num].set_ylim([0, 1.1*peak_value])

    elif num_to_plot == 2:
        letters = ["a)", "b)"]
        fig, ax = plt.subplots(1,2)
        peak_value = -np.inf

        for i, run_name in enumerate(run_names):

            # load results
            save_folder = os.path.join(cwd, "output", run_name)
            disease_one_counts = np.load(os.path.join(save_folder, "disease_one_counts.npy"))
            disease_two_counts = np.load(os.path.join(save_folder, "disease_two_counts.npy"))
            ax[i].plot(disease_one_counts[:,1], label="Disease 1")
            ax[i].plot(disease_two_counts[:,1], label="Disease 2")
            ax[i].grid(True)
            ax[i].legend()
            ax[i].set_xlabel("Day")

            # adjust calculated peak for potential standardization
            peak_value = max(peak_value, np.max(disease_one_counts[:,1]), np.max(disease_two_counts[:,2]))

            if i == 0:
                ax[i].set_ylabel("Infected nodes")               

            # add letterin
            ax[i].text(
                    -0.05, 1.07, letters[i], transform=(
                        ax[i].transAxes + ScaledTranslation(-25/72, -25/72, fig.dpi_scale_trans)),
                    fontsize='large', va='bottom', fontfamily='serif')

            # adjust caption
            if i==0:
                caption += f"{letters[i]} {run_names[i]} and "
            else:
                caption += f"{letters[i]} {run_names[i]}."
        
        if standardize_y_axis:
            for i in range(0, 2):
                ax[i].set_ylim([0, 1.1*peak_value])
    else:
        raise ValueError("Number of plots not implemented")

    print(caption)
    plt.show()

if __name__ == "__main__":
    comparative_outputs(["less_extreme_response", "low_extreme_response"], standardize_y_axis=True)
    #comparative_outputs(["no_recovery", "recovery_0.5", "recovery_1", "recovery_2"], standardize_y_axis=True)