from simulation import *
import matplotlib.pyplot as plt

# run_dict = {"no_prior_res": {"prior_res": [0,0]},
#             "prior_res_0.2": {"prior_res": [0,0.2]},
#             "prior_res_0.4": {"prior_res": [0,0.4]}}

# for run_name in run_dict.keys():
#     run_simulation(run_name, opt_args=run_dict[run_name])

# run_simulation("retrieve_new_setup", opt_args={"rec_weights": [1,1]}, node_folder="default_matrix")
# plot_outputs("retrieve_new_setup")

run_name = "same_low_virulence_long_illness"
#opt_args = {"trans_prob": [0.01, 0.01], "fall_ill": [0.05, 0.1]}
#run_simulation(run_name, opt_args=opt_args, node_folder="default_matrix")
#plot_outputs(run_name)
#plot_outputs("same_virulence")
to_redo = "less_connectivity"
re_run_simulation(to_redo)
plot_outputs(to_redo)

# currently valid
# default_matrix, less_connectivity, less_extreme_response, low_extreme_response