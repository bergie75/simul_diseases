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

if __name__ == "__main__":
    run_names = []
    delete_runs(run_names)