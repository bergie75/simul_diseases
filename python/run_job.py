from simulation import *

run_name = "default_two_week_sine"
opt_args = {"trans_mod": "sine_two_week"}
run_simulation(run_name, opt_args=opt_args, node_folder="default_matrix")
plot_outputs(run_name)

# to_redo = "less_connectivity"
# re_run_simulation(to_redo)
# plot_outputs(to_redo)
