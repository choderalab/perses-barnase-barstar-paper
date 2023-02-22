import os
import pickle
import argparse
from pathlib import Path

from openmm import unit
from perses.annihilation.relative import HybridTopologyFactory

# Read args
parser = argparse.ArgumentParser(description='run equilibration')
parser.add_argument('in_dir', type=str, help='input directory')
parser.add_argument('out_dir', type=str, help='output directory')
args = parser.parse_args()

phases = ['complex', 'apo', 'vacuum']

for phase in phases:
    directory_number = Path(args.in_dir).parts[-2]
    if os.path.exists(os.path.join(args.in_dir, f"{directory_number}_{phase}.pickle")):

        print(f"GENERATING HTF FOR PHASE: {phase}")

        # Load old htf
        with open(os.path.join(args.in_dir, f"{directory_number}_{phase}.pickle"), 'rb') as f:
            htf = pickle.load(f) 

        # Create vanilla htf
        vanilla_htf = HybridTopologyFactory(topology_proposal=htf._topology_proposal,
                         current_positions=htf.old_positions(htf.hybrid_positions),
                         new_positions=htf.new_positions(htf.hybrid_positions),
                         softcore_LJ_v2=False,
                         interpolate_old_and_new_14s=True,
                         impose_virtual_bonds=False
                        )

        # Check the softcore potential
        force_dict = {force.getName(): index for index, force in enumerate(vanilla_htf.hybrid_system.getForces())}
        print("CustomNonbondedForce expression: ", vanilla_htf.hybrid_system.getForce(force_dict['CustomNonbondedForce']).getEnergyFunction())

        # Save htf
        directory_number = Path(args.out_dir).parts[-2]
        with open(os.path.join(args.out_dir, f"{directory_number}_{phase}.pickle"), "wb") as f:
            pickle.dump(vanilla_htf, f)

