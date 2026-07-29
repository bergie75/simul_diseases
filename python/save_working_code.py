    # parameter folder
    cwd = os.getcwd()
    par_folder = os.path.join(cwd, "parameters")
    new_matrix = False

    # other parameters
    num_days = 90
    scale = 10**3
    num_nodes = 40*scale
    total_node_set = set(range(0,num_nodes))
    starting_exposure_frac = np.array([1.0, 1.0])/scale
    delay_day = 0

    # will keep track of all daily statistics
    daily_counts = [[], []]
    affected_nodes = [set(), set()]
    global_blocking_events = []
    count_self_blocks = True
    used_recovered_node = np.zeros((2,num_days))
    used_recovered_degree = np.zeros((2,num_days))


    # delayed entry of second disease options
    if delay_day > 0:
        starting_exposure_frac[-1] = 0
    late_frac = 1.0/scale

    num_diseases = len(starting_exposure_frac)
    trans_prob = [0.04, 0.01]
    fall_ill = [0.2,0.1]  # prob of becoming fully infected
    rec_prob = [0.5/num_days,0.5/num_days]
    dis_weight = [1,1]
    prior_res = [0,0]

    # adjacency matrix construction and/or saving
    if new_matrix:
        #G = newman_watts_strogatz_graph(num_nodes, int(num_nodes/100), 0)
        #G = erdos_renyi_graph(num_nodes, 0.7/scale)
        # G=barabasi_albert_graph(num_nodes, int(0.3/scale*num_nodes))
        # adjacency_matrix = nx.to_scipy_sparse_array(G)
        adjacency_matrix = construct_adj_mat(num_nodes, 100, 0.5/scale, 30)
        np.save(os.path.join(par_folder, "adj_mat"), adjacency_matrix)
        
    else:
        adjacency_matrix = np.load(os.path.join(par_folder, "adj_mat.npy"))

    # calculate degree list
    degree_values = np.sum(adjacency_matrix, axis=0)
    print(f"Maximum degree: {max(degree_values)}")
    print(f"Mean degree value: {np.mean(degree_values)}")

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
        pos = nx.spring_layout(G, seed=42335487)
        nx.draw_networkx_nodes(G, pos, nodelist=list(affected_nodes[0]), node_color="blue")
        nx.draw_networkx_nodes(G, pos, nodelist=list(affected_nodes[1]), node_color="orange")
        remaining_nodes = total_node_set.difference(affected_nodes[0].union(affected_nodes[1]))
        nx.draw_networkx_nodes(G, pos, nodelist=list(remaining_nodes), node_color="gray")
        nx.draw_networkx_edges(G, pos, width=0.3, alpha=0.5)

    # create subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)
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
            if status[1-disease_index] == "E":
                exposure_blocks[disease_index] += 1
            elif status[1-disease_index] == "I":
                infection_blocks[disease_index] += 1

        if not block_issued and status[1-disease_index] == "R":
            used_recovered_node[disease_index, day] += 1
            used_recovered_degree[disease_index, day] += degree_values[index]

    ax1.plot(disease_one_counts[:,1], label="Disease 0 infected")
    ax1.plot(disease_two_counts[:,1], label="Disease 1 infected")
    # ax2.plot(disease_one_counts[:,2], label="Disease 0 infected")
    # ax2.plot(disease_two_counts[:,2], label="Disease 1 infected")
    ax2.plot(disease_one_counts[:,3], label="Disease 0 recovered")
    ax2.plot(disease_two_counts[:,3], label="Disease 1 recovered")
    ax3.plot(used_recovered_node[0,:], label="Disease 0 recovered nodes used")
    ax3.plot(used_recovered_node[1,:], label="Disease 1 recovered nodes used")
    ax4.plot(used_recovered_degree[0,:], label="Disease 0 recovered degree used")
    ax4.plot(used_recovered_degree[1,:], label="Disease 1 recovered degree used")   
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