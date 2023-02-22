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
parser.add_argument('replica', type=int, help='replica')

args = parser.parse_args()

def get_residue_distances(phase, out_dir, sub_dir, replica, mutating_residue_index, n_iterations_start=None, n_iterations_end=None, barstar_first=True):
    
    """
    Get distances between residue pairs (distance values are in nm)

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
    mutating_residue_index : int
        mutating residue index
    n_iterations_start : int, default None
        iteration to start the trajectory
    n_iterations_end : int, default None
        iteration to end the trajectory
    barstar_first : booleana, default True
        whether barstar is first in this topology (i.e., the mutation is in barstar)

    Returns
    -------
    data : list of lists
        each sub list contains [replica, iteration, residue pair distance 1, residue pair distance 2, ...]
    column_names : list of str
        names for each corresponding value in data

    """
    
    from mdtraj import compute_contacts
    
    total_iterations = n_iterations_end - n_iterations_start
    distances_all = []
    
    ## OLD POSITIONS ##
    position_type = 'old'
    
    # Load pdb and traj
    pdb_old = md.load(os.path.join(out_dir, f"{sub_dir}_{phase}_{position_type}_replica_{replica}_{n_iterations_start}_{n_iterations_end}_no_imaging.pdb"))
    traj_old = md.load_dcd(os.path.join(out_dir, f"{sub_dir}_{phase}_{position_type}_replica_{replica}_{n_iterations_start}_{n_iterations_end}_no_imaging.dcd"), top=pdb_old)

    # Assemble interface residue pair list
    barstar_chain = 0 if barstar_first else 2
    barnase_chain = 2 if barstar_first else 0
    d_resSeq_to_index_barstar = {res.resSeq: res.index for res in pdb_old.topology.residues if res.chain.index == barstar_chain}
    d_resSeq_to_index_barnase = {res.resSeq: res.index for res in pdb_old.topology.residues if res.chain.index == barnase_chain}
    barnase_indices = [d_resSeq_to_index_barnase[res] for res in barnase_residues]
    barstar_indices = [d_resSeq_to_index_barstar[res] for res in barstar_residues]
    contacts = [pair for pair in itertools.combinations(barnase_indices + barstar_indices, 2)]

    # Retrieve distances
    distances_old, residue_pairs_old = compute_contacts(traj_old, contacts, scheme='closest-heavy')
    
    ## NEW POSITIONS ##
    position_type = 'new'
    
    # Load pdb and traj
    pdb_new = md.load(os.path.join(out_dir, f"{sub_dir}_{phase}_{position_type}_replica_{replica}_{n_iterations_start}_{n_iterations_end}_no_imaging.pdb"))
    traj_new = md.load_dcd(os.path.join(out_dir, f"{sub_dir}_{phase}_{position_type}_replica_{replica}_{n_iterations_start}_{n_iterations_end}_no_imaging.dcd"), top=pdb_new)

    # Retrieve residue pairs involving mutating residue
    contacts_mutating = [pair for pair in contacts if pair[0] == mutating_residue_index or pair[1] == mutating_residue_index]
    
    # Retrieve distances
    distances_new, residue_pairs_new = compute_contacts(traj_new, contacts_mutating, scheme='closest-heavy')
    
    # Concatenate old and new distances
    distances_all = np.concatenate([distances_old, distances_new], axis=1)
    
    # Create rows for dataframe
    data = []
    for iteration in range(total_iterations):
        row = [replica, iteration]
        for i, pair in enumerate(contacts + contacts_mutating):
            row.append(distances_all[iteration][i])
        data.append(row)
    
    # Create column names
    column_names = ["replica", "iteration"]
    d_chain_name = {barstar_chain: "barstar", barnase_chain: "barnase"}
    d_index_to_name = {res.index: f"{d_chain_name[res.chain.index]}-{res.resSeq}" for res in pdb_old.topology.residues if res.chain.index in [0, 2]}
    for pair in contacts:
        column_names.append(f"{d_index_to_name[pair[0]]} | {d_index_to_name[pair[1]]} old")
    for pair in contacts_mutating:
        column_names.append(f"{d_index_to_name[pair[0]]} | {d_index_to_name[pair[1]]} new")
    
    return data, column_names

## Define interface residues
# Load crystal structure in 
barstar = md.load_pdb("/data/chodera/zhangi/perses_benchmark/perses-barnase-barstar-paper/input_files/1brs_barstar_renumbered.pdb")
barnase = md.load_pdb("/data/chodera/zhangi/perses_benchmark/perses-barnase-barstar-paper/input_files/1brs_barnase_renumbered.pdb")

bnbs_topology = barnase.topology.join(barstar.topology)
bnbs_positions = np.concatenate((barnase.xyz, barstar.xyz), axis=1)

bnbs = md.Trajectory(xyz=bnbs_positions, topology=bnbs_topology)

# Get atom indices for barnase and barstar
barnase_atoms = bnbs.topology.select("chainid == 0 and protein")
barstar_atoms = bnbs.topology.select("chainid == 2 and protein")

# Set radius (nm)
radius = 0.4

# Retrieve barnase residues within 4 angstroms of the barstar
neighbors_barnase_atoms = md.compute_neighbors(bnbs, radius, barstar_atoms, haystack_indices=barnase_atoms)[0]
neighbors_barnase_residues = set([atom.residue.resSeq for atom in bnbs.topology.atoms if atom.index in neighbors_barnase_atoms])

# Retrieve barstar residues within 4 angstroms of the barnase
neighbors_barstar_atoms = md.compute_neighbors(bnbs, radius, barnase_atoms, haystack_indices=barstar_atoms)[0]
neighbors_barstar_residues = set([atom.residue.resSeq for atom in bnbs.topology.atoms if atom.index in neighbors_barstar_atoms])

# Print barnase interface residues and create atom set dict
barnase_residues = []
for res in bnbs.topology.residues:
    if res.chain.index == 0 and res.resSeq in neighbors_barnase_residues:
        print(res)
        barnase_residues.append(res.resSeq)

# Print barstar interface residues and create atom set dict
barstar_residues = []
for res in bnbs.topology.residues:
    if res.chain.index == 2 and (res.resSeq in neighbors_barstar_residues or res.resSeq == 80):
        print(res)
        barstar_residues.append(res.resSeq)
    
# Set parameters
phase = args.phase
out_dir = args.out_dir
sub_dir = args.sub_dir
main_dir = int(out_dir.split('/')[-3])
n_iterations_start = 0
if main_dir == 45:
    n_replicas = 24 if sub_dir in [9, 11] else 36
    barstar_first = False if sub_dir in [10] else True
    n_iterations_end = 10001 if sub_dir in [9, 10] else 1001
    # Define mapping between sub_dir and mutating residue index
    d_mutating_residue_index = {9: 42, 10: 85, 11: 44}
else:
    n_replicas = 24 if sub_dir in [0, 1, 2, 3, 4, 14, 15, 16, 17, 18] else 36
    barstar_first = False if sub_dir in [9, 10, 11, 12, 13, 23, 24, 25, 26, 27] else True
    n_iterations_end = 10001 if sub_dir in [12, 17] else 1001
    # Define mapping between sub_dir and mutating residue index
    d_mutating_residue_index = {0: 29, 1: 29, 2: 38, 3: 42, 4: 44, 5: 76, 6: 80, 7: 35, 8: 39,
                            9: 25, 10: 57, 11: 81, 12: 85, 13: 100,
                            14: 29, 15: 29, 16: 38, 17: 42, 18: 44, 19: 76, 20: 80, 21: 35, 22: 39,
                            23: 25, 24: 57, 25: 81, 26: 85, 27: 100
                           }
mutating_residue_index = d_mutating_residue_index[sub_dir]

print(f"out_dir: {out_dir}")
print(f"mutating residue index: {mutating_residue_index}")
print(f"n_replica: {n_replicas}")
print(f"barstar_first? {barstar_first}")

# Get distances and save dataframe
for replica in tqdm_notebook(range(n_replicas)):
    data, column_names = get_residue_distances(phase, out_dir, sub_dir, replica, mutating_residue_index, n_iterations_start, n_iterations_end, barstar_first=barstar_first)
    df = pd.DataFrame(data, columns=column_names)
    with open(os.path.join(out_dir, f"dataframe_residue_distances_replica{replica}.pickle"), "wb") as f:
        pickle.dump(df, f)
    
