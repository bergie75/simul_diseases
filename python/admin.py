import os
import shutil

def delete_runs(run_names):
    cwd = os.getcwd()
    for run_name in run_names:
        par_folder = os.path.join(cwd, "parameters", run_name)
        output_folder = os.path.join(cwd, "output", run_name)

        if os.path.isdir(par_folder):
            shutil.rmtree(par_folder)

        if os.path.isdir(output_folder):
            shutil.rmtree(output_folder)

def rename_runs(old_run_names, new_run_names):
    cwd = os.getcwd()
    for i,run_name in enumerate(old_run_names):
        old_par_folder = os.path.join(cwd, "parameters", run_name)
        old_output_folder = os.path.join(cwd, "output", run_name)

        if os.path.isdir(old_par_folder):
            new_par_folder  = os.path.join(cwd, "parameters", new_run_names[i])
            os.rename(old_par_folder, new_par_folder)

        if os.path.isdir(old_output_folder):
            new_output_folder = os.path.join(cwd, "output", new_run_names[i])
            os.rename(old_output_folder, new_output_folder)

if __name__ == "__main__":
    run_names = []
    delete_runs(run_names)

    old_run_names = ["default_two_week_trans_cycle"]
    new_run_names = ["default_two_week_cosine"]
    rename_runs(old_run_names, new_run_names)