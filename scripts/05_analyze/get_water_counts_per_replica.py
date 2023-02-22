import os
import sys
import numpy as np
import mdtraj as md
import pickle
from simtk.openmm import unit, app
from tqdm import tqdm_notebook
from openeye import oechem
import pandas as pd
import itertools

import argparse

# Load arguments
parser = argparse.ArgumentParser(description='run repex')
parser.add_argument('out_dir', type=str, help='path to input/output dir')
parser.add_argument('sub_dir', type=int, help='sub directory')
parser.add_argument('phase', type=str, help='phase')

args = parser.parse_args()

def get_nearby_waters(phase, out_dir, sub_dir, replica, mutating_residue_id, mutating_residue_names, radius=0.5, n_iterations_start=None, n_iterations_end=None):
    
    """
    Get number of waters within 5 angstroms of the mutating residue (distance values are in nm)

    Parameters
    ----------
    phase : str
        phase (e.g., 'complex')
    out_dir : str
        path to directory containing data to analyze
    sub_dir : str
        path to sub-directory containing data to analyze
    replica : int
        replica id
    mutating_residue_id : int
        mutating residue id
    mutating residue names : list of str
        the three letter codes of the WT and mutant residues
    radius : float
        radius of neighborhood (in nanometers)
    n_iterations_start : int, default None
        iteration to start the trajectory
    n_iterations_end : int, default None
        iteration to end the trajectory

    Returns
    -------
    data : list of lists
        each sub list contains [replica, iteration, neighboring waters for old topology, neighboring waters for new topology]
    column_names : list of str
        names for each corresponding value in data
    
    """
    
    from mdtraj import compute_neighbors
    
    total_iterations = n_iterations_end - n_iterations_start
    
    neighbor_counts_all = []
    position_types = ['old', 'new']
    for position_type, mutating_residue_name in zip(position_types, mutating_residue_names): 
    
        # Load pdb and traj
        pdb = md.load(os.path.join(out_dir, f"{sub_dir}_{phase}_{position_type}_replica_{replica}_{n_iterations_start}_{n_iterations_end}_no_imaging.pdb"))
        traj = md.load_dcd(os.path.join(out_dir, f"{sub_dir}_{phase}_{position_type}_replica_{replica}_{n_iterations_start}_{n_iterations_end}_no_imaging.dcd"), top=pdb)

        # Retrieve heavy atoms in mutating residue
        target_atoms = traj.topology.select(f"resSeq == {mutating_residue_id} and resname == {mutating_residue_name} and element != H")
        for atom in traj.topology.atoms:
            if atom.index in target_atoms:
                    print(atom, atom.residue.resSeq, atom.index)

        # Get waters near target atoms
        water_atoms = traj.topology.select("water and element == O")
        neighbors = md.compute_neighbors(traj, radius, target_atoms, haystack_indices=water_atoms)
        neighbor_counts = [len(neighbors_per_frame) for neighbors_per_frame in neighbors]
        neighbor_counts_all.append(neighbor_counts)
    
    # Create rows for dataframe
    data = []
    for iteration in range(total_iterations):
        row = [replica, iteration]
        for index, position_type in enumerate(neighbor_counts_all):
            row.append(neighbor_counts_all[index][iteration])
        data.append(row)
    
    # Create column names
    column_names = ["replica", "iteration", "waters nearby mutating res (radius = 0.5 nm) old", "waters nearby mutating res (radius = 0.5 nm) new"]
    
    return data, column_names

# Set parameters
phase = args.phase
out_dir = args.out_dir
sub_dir = args.sub_dir
main_dir = int(out_dir.split('/')[-3])
n_iterations_start = 0
if main_dir == 45:
    n_replicas = 24 if sub_dir in [9, 11] else 36
    n_iterations_end = 10001 if sub_dir in [9, 10] else 1001
    # Define mapping between sub_dir and mutating residue id
    d_mutating_resnames = {9: ["ALA", "THR"], 10: ["ARG", "ALA"], 11: ["TRP", "PHE"]}
    # Define mapping between sub_dir and mutating residue names
    d_mutating_resid = {9: 42, 10: 87, 11: 44}
else:
    n_replicas = 24 if sub_dir in [0, 1, 2, 3, 4, 14, 15, 16, 17, 18] else 36
    n_iterations_end = 10001 if sub_dir in [12, 17] else 1001
    # Define mapping between sub_dir and mutating residue id
    d_mutating_resid = {0: 29, 1: 29, 2: 38, 3: 42, 4: 44, 5: 76, 6: 80, 7: 35, 8: 39,
                            9: 27, 10: 59, 11: 83, 12: 87, 13: 102,
                            14: 29, 15: 29, 16: 38, 17: 42, 18: 44, 19: 76, 20: 80, 21: 35, 22: 39,
                            23: 27, 24: 59, 25: 83, 26: 87, 27: 102
                           }
    # Define mapping between sub_dir and mutating residue names
    d_mutating_resnames = {0: ["TYR", "PHE"], 1: ["TYR", "ALA"], 2: ["TRP", "PHE"], 3: ["THR", "ALA"], 4: ["TRP", "PHE"],
                           5: ["GLU", "ALA"], 6: ["GLU", "ALA"], 7: ["ASP", "ALA"], 8: ["ASP", "ALA"],
                           9: ["LYS", "ALA"], 10: ["ARG", "ALA"], 11: ["ARG", "GLN"], 12: ["ARG", "ALA"], 13: ["HIS", "ALA"],
                           14: ["PHE", "TYR"], 15: ["ALA", "TYR"], 16: ["PHE", "TRP"], 17: ["ALA", "THR"],  18: ["PHE", "TRP"], 
                           19: ["ALA", "GLU"], 20: ["ALA", "GLU"], 21: ["ALA", "ASP"], 22: ["ALA", "ASP"],
                           23: ["ALA", "LYS"], 24: ["ALA", "ARG"], 25: ["GLN", "ARG"], 26: ["ALA", "ARG"], 
                           27: ["ALA", "HIS"],
                           }
mutating_residue_id = d_mutating_resid[sub_dir]
mutating_residue_names = d_mutating_resnames[sub_dir]

print(f"out_dir: {out_dir}")
print(f"mutating residue id: {mutating_residue_id}")
print(f"mutating residue names: {mutating_residue_names}")
print(f"n_replicas: {n_replicas}")

# Get water counts and save dataframe
data_all_replicas = []
for replica in tqdm_notebook(range(n_replicas)):
    data, column_names = get_nearby_waters(phase, out_dir, sub_dir, replica, mutating_residue_id, mutating_residue_names, n_iterations_start=n_iterations_start, n_iterations_end=n_iterations_end)
    data_all_replicas.append(data)

data_all_replicas_final = []
for data in data_all_replicas:
    for row in data:
        data_all_replicas_final.append(row)

df = pd.DataFrame(data_all_replicas_final, columns=column_names)

with open(os.path.join(out_dir, "dataframe_waters_per_replica.pickle"), "wb") as f:
    pickle.dump(df, f) 
