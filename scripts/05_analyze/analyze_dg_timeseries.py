import os
import sys
import pickle
import argparse
from tqdm import tqdm
from analysis_tools import DataAnalyzer

# Load arguments
parser = argparse.ArgumentParser(description='analyze dg')
parser.add_argument('main_dir', type=int, help='main dir')
parser.add_argument('sub_dir', type=int, help='sub dir')
parser.add_argument('replicate', type=int, help='replicate number')
parser.add_argument('phase', type=str, help='phase of the simulation to use in storage filename')
parser.add_argument('starting_iteration', type=int, help='starting iteration to analyze' )
parser.add_argument('total_iterations', type=int, help='total number of iterations to analyze')
parser.add_argument('--initialize_with_bar', action='store_true', help='whether to compute free energy with the BAR estimate as the initial estimate')
args = parser.parse_args()

# Set parameters
main_dir = args.main_dir
replicate = args.replicate
phase = args.phase
sub_dir = args.sub_dir
path = f"/data/chodera/zhangi/perses_benchmark/repex/perses-bnbs-paper-fourth-attempt/{main_dir}/{sub_dir}/replicate_{replicate}/{sub_dir}_{phase}.nc"
analysis_interval = 1000 # iterations
discard_fraction = 0.1
subsample_rate = 5

# Set MBAR initialization
initialize = 'zeros' if not args.initialize_with_bar else 'BAR'
print(f"Free energy estimates will be initalized with {initialize}")

# Create analyzer
data_analyzer = DataAnalyzer(path) 

# Retrieve free energies
results_all = []
metadata_all = []
for max_n_iterations in range(args.starting_iteration, args.total_iterations + 1, analysis_interval):
    print(max_n_iterations, file=sys.stderr)
    results, metadata = data_analyzer.get_free_energy(bootstrap=True, max_n_iterations=max_n_iterations, n_equilibration_iterations=round(max_n_iterations*discard_fraction), subsample_rate=subsample_rate, initialize=initialize)
    results_all.append(results)
    metadata_all.append(metadata)

# Save data as pickles
starting_ns = int(args.starting_iteration/1000) # 1000 iterations = 1 ns
total_ns = int(args.total_iterations/1000) # 1000 iterations = 1 ns
with open(os.path.join(os.path.dirname(path), f"{phase}_free_energy_timeseries_{starting_ns}_{total_ns}ns.pickle"), "wb") as f:
    pickle.dump(results_all, f)
with open(os.path.join(os.path.dirname(path), f"{phase}_free_energy_timeseries_{starting_ns}_{total_ns}ns_metadata.pickle"), "wb") as f:
    pickle.dump(metadata_all, f)

