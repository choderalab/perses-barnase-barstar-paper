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

def get_dihedral_angles(phase, out_dir, sub_dir, replica, n_iterations_start=None, n_iterations_end=None, barstar_first=True):
    
    """
    Angle values are in radians: https://github.com/mdtraj/mdtraj/blob/1a757274f6f9acd6c9c847c72f5857e9c1046d17/mdtraj/geometry/dihedral.py#L105

    Parameters
    ----------
    phase : str
        phase of simulation (e.g., 'complex')
    out_dir : str
        path to directory containing data to analyze
    sub_dir : str
        path to sub-directory containing data to analyze
    replica : int, default None
        replica id
    n_iterations_start : int, default None
        iteration to start the trajectory
    n_iterations_end : int, default None
        iteration to end the trajectory
    barstar_first : boolean, default True
        whether barstar is first in this topology (i.e., the mutation is in barstar)
    
    Returns
    -------
    data : list of lists
        each sublist contains [replica, res, iteration, old phi angle value, old psi angle value,
            old chi1 angle value, old chi2 angle value, old chi3 angle value, old chi4 angle value,
            new phi angle value, new psi angle value, new chi1 angle value, new chi2 angle value,
            new chi3 angle value, new chi4 angle value]

    """
    
    from mdtraj import compute_phi, compute_psi, compute_chi1, compute_chi2, compute_chi3, compute_chi4
    
    def retrieve_angles_of_interest(indices, angles, is_phi_psi=False):
        """
        From mdtraj output, isolate angles for residues of interest

        Parameters
        ----------
        indices : list of lists
            the atom indices for each dihedral computed
        angles : list of lists
            the angle values corresponding to each of the dihedral angles for each frame
        is_phi_psi : boolean, default False
            whether the dihedral angle is phi or psi

        Returns
        -------
        d_angles_of_interest : dict
            angle values for each residue of interest, across all frames
            key : residue name (format: {barstar|barnase}-{residue id}), example: 'barstar-42' 
            value : angle values across all frames

        """
        
        intersection_size = 3 if is_phi_psi else 4 # for phi/psi angles, 3 of the atoms will be from the residue of interest and 1 will be from neighboring residue
        
        d_angles_of_interest = {}
        for res in d_atom_sets:
            atom_set = set(d_atom_sets[res])
            for index, dihedral_indices in enumerate(indices):
                dihedral_indices = set(dihedral_indices)
                if len(dihedral_indices.intersection(atom_set)) == intersection_size:
                    d_angles_of_interest[res] = angles[:, index]
        
        return d_angles_of_interest
    
    total_iterations = n_iterations_end - n_iterations_start
    dicts = []
    for position_type in ['old', 'new']:
    
        # Load pdb and traj
        pdb = md.load(os.path.join(out_dir, f"{sub_dir}_{phase}_{position_type}_replica_{replica}_{n_iterations_start}_{n_iterations_end}_no_imaging.pdb"))
        traj = md.load_dcd(os.path.join(out_dir, f"{sub_dir}_{phase}_{position_type}_replica_{replica}_{n_iterations_start}_{n_iterations_end}_no_imaging.dcd"), top=pdb)

        # Create atom set dict for interface residues
        barstar_chain = 0 if barstar_first else 2
        barnase_chain = 2 if barstar_first else 0
        d_atom_sets = {}
        for res in pdb.topology.residues:
            k = res.resSeq
            if (k in barnase_residues and res.chain.index == barnase_chain):
                d_atom_sets["barnase-" + str(k)] = [atom.index for atom in res.atoms]
            if(k in barstar_residues and res.chain.index == barstar_chain):
                d_atom_sets["barstar-" + str(k)] = [atom.index for atom in res.atoms]

        # Retrieve angles
        for angle_index, compute_angle in enumerate([compute_phi, compute_psi, compute_chi1, 
                                           compute_chi2, compute_chi3, compute_chi4]):
            
            # Get angle indices and angle values for each iteration
            indices, angles = compute_angle(traj)
            
            # Retrieve angles of interest
            if angle_index < 2:
                d_angles_of_interest = retrieve_angles_of_interest(indices, angles, is_phi_psi=True)
            else:
                d_angles_of_interest = retrieve_angles_of_interest(indices, angles)
            
            dicts.append(d_angles_of_interest)

    # Create rows for dataframe
    data = []
    d_phi_angles_of_interest = dicts[0]
    for res in d_phi_angles_of_interest:
        for iteration in range(total_iterations):
            row = [replica, res, iteration]
            for d in dicts:
                if res in d:
                    row.append(d[res][iteration])
                else:
                    row.append(np.nan)
            data.append(row)
                
    return data


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
else:
    n_replicas = 24 if sub_dir in [0, 1, 2, 3, 4, 14, 15, 16, 17, 18] else 36
    barstar_first = False if sub_dir in [9, 10, 11, 12, 13, 23, 24, 25, 26, 27] else True  
    n_iterations_end = 10001 if sub_dir in [12, 17] else 1001

print(f"out_dir: {out_dir}")
print(f"main_dir: {main_dir}")
print(f"n_replicas: {n_replicas}")
print(f"barstar_first? {barstar_first}")

# Get dihedrals and save dataframe
for replica in range(n_replicas):
    data = get_dihedral_angles(phase, out_dir, sub_dir, replica, n_iterations_start, n_iterations_end, barstar_first=barstar_first)
    df = pd.DataFrame(data, columns=["replica", "residue", "iteration", 
                                 "phi_old", "psi_old", "chi1_old", "chi2_old", "chi3_old", "chi4_old", 
                                 "phi_new", "psi_new", "chi1_new", "chi2_new", "chi3_new", "chi4_new"])
    with open(os.path.join(out_dir, f"dataframe_dihedral_angles_per_replica{replica}.pickle"), "wb") as f:
        pickle.dump(df, f)
