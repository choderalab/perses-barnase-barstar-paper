import os
import pickle
import argparse
from pathlib import Path

from openmm import unit
from perses.annihilation.relative import RESTCapableHybridTopologyFactory

# Read args
parser = argparse.ArgumentParser(description='run equilibration')
parser.add_argument('in_dir', type=str, help='input directory from which to copy the htfs')
parser.add_argument('out_dir', type=str, help='output directory')
parser.add_argument('--w_lifting', type=float, default=0.3, help='w_lifting to use when generating modified htf')
parser.add_argument('--rest_radius', type=float, default=0.3, help='rest_radius to use when generating modified htf')
args = parser.parse_args()

phases = ['complex', 'apo', 'vacuum']

for phase in phases:
    directory_number = Path(args.in_dir).parts[-2]
    out_directory_number = Path(args.out_dir).parts[-2]
    if os.path.exists(os.path.join(args.in_dir, f"{directory_number}_{phase}.pickle")):

        # Load old htf
        with open(os.path.join(args.in_dir, f"{directory_number}_{phase}.pickle"), 'rb') as f:
            htf = pickle.load(f) 

        print(f"GENERATING HTF WITH W_LIFTING: {args.w_lifting} AND REST_RADIUS: {args.rest_radius} FOR PHASE: {phase}")

        # Create modified htf
        modified_htf = RESTCapableHybridTopologyFactory(topology_proposal=htf._topology_proposal,
                         current_positions=htf.old_positions(htf.hybrid_positions),
                         new_positions=htf.new_positions(htf.hybrid_positions),
                         rest_radius=args.rest_radius * unit.nanometers,
                         w_lifting=args.w_lifting * unit.nanometers
                        )

        # Save htf
        with open(os.path.join(args.out_dir, f"{out_directory_number}_{phase}.pickle"), "wb") as f:
            pickle.dump(modified_htf, f)

