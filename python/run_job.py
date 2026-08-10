from simulation import *
import matplotlib.pyplot as plt

# run_dict = {"no_prior_res": {"prior_res": [0,0]},
#             "prior_res_0.2": {"prior_res": [0,0.2]},
#             "prior_res_0.4": {"prior_res": [0,0.4]}}

# for run_name in run_dict.keys():
#     run_simulation(run_name, opt_args=run_dict[run_name])

# run_simulation("retrieve_new_setup", opt_args={"rec_weights": [1,1]}, node_folder="default_matrix")
# plot_outputs("retrieve_new_setup")

run_name = "less_extreme_response"
opt_args = {"large_resp":[0.5,0.5]}
run_simulation(run_name, opt_args=opt_args, node_folder="default_matrix")
plot_outputs(run_name)
# to_redo = "less_extreme_response"
# re_run_simulation(to_redo)
# plot_outputs(to_redo)

# currently valid
# default_matrix, less_connectivity, less_extreme_response, low_extreme_response