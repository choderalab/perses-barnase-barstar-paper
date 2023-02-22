import os
import sys
import numpy as np
import mdtraj as md
import pickle
from simtk.openmm import unit, app
from tqdm import tqdm_notebook
from openeye import oechem

import argparse

# Load arguments
parser = argparse.ArgumentParser(description='run repex')
parser.add_argument('out_dir', type=str, help='path to input/output dir')
parser.add_argument('sub_dir', type=str, help='phase of the simulation to use in storage filename')
parser.add_argument('phase', type=str, help='number of states')
parser.add_argument('replica', type=int, help='number of iterations to run')
args = parser.parse_args()

def make_traj(phase, out_dir, sub_dir, htf, state=None, replica_id=None, n_iterations_start=None, n_iterations_end=None):
    """
    Retrieve pdbs/dcds of the old and new positions for a given thermodynamic state index. 
    
    Adapted from Hannah: https://github.com/hannahbrucemacdonald/endstate_pdbs/blob/master/scripts/input_for_pol_calc.py

    Parameters
    ----------
    phase : str
        phase of simulation (e.g., 'complex')
    out_dir : str
        path to directory containing data to analyze
    sub_dir : str
        path to sub-directory containing data to analyze
    htf : perses.annihilation.relative.RESTCapableHybridTopologyFactory
        hybrid topology factory containing hybrid topology, positions, system
    state : int, default None
        state index
    replica_id : int, default None
        replica id
    n_iterations_start : int, default None
        iteration to start the trajectory
    n_iterations_end : int, default None
        iteration to end the trajectory

    """
    
    # Create mdtraj topologies that correspond to the positions we want to create trajectories for
    new_top = md.Topology.from_openmm(htf._topology_proposal.new_topology)
    old_top = md.Topology.from_openmm(htf._topology_proposal.old_topology)
    
    # Load nc files
    from perses.analysis.utils import open_netcdf
    nc = open_netcdf(os.path.join(out_dir, f"{sub_dir}_{phase}.nc"))
    nc_checkpoint = open_netcdf(os.path.join(out_dir, f"{sub_dir}_{phase}_checkpoint.nc"))
    checkpoint_interval = nc_checkpoint.CheckpointInterval
    all_positions = nc_checkpoint.variables['positions']
    box_vectors = np.array(nc_checkpoint['box_vectors'])
    n_iter, n_replicas, n_atoms, _ = np.shape(all_positions)
    n_iterations_start = 0 if n_iterations_start is None else int(n_iterations_start / checkpoint_interval)
    n_iterations_end = n_iter if n_iterations_end is None else int(n_iterations_end / checkpoint_interval)
    n_iterations = n_iterations_end - n_iterations_start
    print(n_iterations_start, n_iterations_end, file=sys.stderr)
    
    # Retrieve positions 
    all_pos_new = np.zeros(shape=(n_iterations, new_top.n_atoms, 3))
    all_pos_old = np.zeros(shape=(n_iterations, old_top.n_atoms, 3))
    all_box_vectors = np.zeros(shape=(n_iterations, 3, 3), dtype=np.float32)
    for iteration in tqdm_notebook(range(n_iterations_start, n_iterations_end)):
        print(iteration, file=sys.stderr)
        if state is not None:
            replica_id = np.where(nc.variables['states'][iteration*checkpoint_interval] == state)[0][0]
        pos = all_positions[iteration,replica_id,:,:] * unit.nanometers
        all_pos_new[iteration-n_iterations_start] = htf.new_positions(pos).value_in_unit_system(unit.md_unit_system)
        all_pos_old[iteration-n_iterations_start] = htf.old_positions(pos).value_in_unit_system(unit.md_unit_system)
        all_box_vectors[iteration-n_iterations_start] = box_vectors[iteration,replica_id,:,:]
    
    # Create trajectories
    traj_old = md.Trajectory(all_pos_old, old_top)
    traj_new = md.Trajectory(all_pos_new, new_top)
    
    # Set unit cell vectors in traj 
    traj_old.unitcell_vectors = all_box_vectors
    traj_new.unitcell_vectors = all_box_vectors
    
    # Save old trajectory
    print("saving old traj", file=sys.stderr)
    traj_old.save(os.path.join(out_dir, f"{sub_dir}_{phase}_old_replica_{replica}_{n_iterations_start}_{n_iterations_end}_no_imaging.dcd"))
    traj_old[0].save(os.path.join(out_dir, f"{sub_dir}_{phase}_old_replica_{replica}_{n_iterations_start}_{n_iterations_end}_no_imaging.pdb"))
    
    # Save new trajectory
    print("saving new traj", file=sys.stderr)
    traj_new.save(os.path.join(out_dir, f"{sub_dir}_{phase}_new_replica_{replica}_{n_iterations_start}_{n_iterations_end}_no_imaging.dcd"))
    traj_new[0].save(os.path.join(out_dir, f"{sub_dir}_{phase}_new_replica_{replica}_{n_iterations_start}_{n_iterations_end}_no_imaging.pdb"))

# Set args
out_dir = args.out_dir
sub_dir = args.sub_dir
phase = args.phase
replica = args.replica

# Load htf
with open(os.path.join(out_dir, f"{sub_dir}_{phase}.pickle"), "rb") as f:
    htf = pickle.load(f)

# Generate old and new trajs
make_traj(phase, out_dir, sub_dir, htf, replica_id=replica)


