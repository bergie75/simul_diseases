from simulation import *
import matplotlib.pyplot as plt

# run_dict = {"no_prior_res": {"prior_res": [0,0]},
#             "prior_res_0.2": {"prior_res": [0,0.2]},
#             "prior_res_0.4": {"prior_res": [0,0.4]}}

# for run_name in run_dict.keys():
#     run_simulation(run_name, opt_args=run_dict[run_name])

# run_simulation("retrieve_new_setup", opt_args={"rec_weights": [1,1]}, node_folder="default_matrix")
# plot_outputs("retrieve_new_setup")

opt_args = {"rec_weights": [1,1], "scaled_p_ER": 0.05, "m_BA": 10}
run_simulation("low_connectivity_less_recovery", opt_args=opt_args, node_folder="low_connectivity")
plot_outputs("low_connectivity_less_recovery")