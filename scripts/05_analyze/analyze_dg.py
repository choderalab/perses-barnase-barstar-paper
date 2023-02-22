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
parser.add_argument('total_iterations', type=int, help='total number of iterations to analyze')
parser.add_argument('--initialize_with_bar', action='store_true', help='whether to compute free energy with the BAR estimate as the initial estimate')
parser.add_argument('--bootstrap', action='store_true', help='whether to bootstrap the free energy (and use a less stringent tolerance)')
args = parser.parse_args()

# Set path
main_dir = args.main_dir
replicate = args.replicate
phase = args.phase
sub_dir = args.sub_dir
path = f"/data/chodera/zhangi/perses_benchmark/repex/perses-bnbs-paper-fourth-attempt/{main_dir}/{sub_dir}/replicate_{replicate}/{sub_dir}_{phase}.nc"

# Set MBAR initialization
initialize = 'zeros' if not args.initialize_with_bar else 'BAR'
print(f"Free energy estimates will be initalized with {initialize}")

# Set bootstrap
bootstrap = args.bootstrap
print(f"Free energy estimates bootstrapped? {bootstrap}")

# Create analyzer
data_analyzer = DataAnalyzer(path) 

# Retrieve free energies
results, metadata = data_analyzer.get_free_energy(bootstrap=bootstrap, max_n_iterations=args.total_iterations, initialize=initialize)

# Save data as pickles
total_ns = int(args.total_iterations/1000) # 1000 iterations = 1 ns
bootstrap_name = '_bootstrapped' if bootstrap else ''
with open(os.path.join(os.path.dirname(path), f"{phase}_free_energy{bootstrap_name}_{total_ns}ns.pickle"), "wb") as f:
    pickle.dump(results, f)
with open(os.path.join(os.path.dirname(path), f"{phase}_free_energy{bootstrap_name}_{total_ns}ns_metadata.pickle"), "wb") as f:
    pickle.dump(metadata, f)

