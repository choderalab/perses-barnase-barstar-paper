from simtk.openmm import unit

from openmmtools.constants import kB
from openmmtools import cache
from openmmtools.multistate import MultiStateReporter
from openmmtools import states

from perses.analysis.utils import open_netcdf
from perses.dispersed.utils import configure_platform

import numpy as np
import os
import copy
import sys

import argparse

platform = configure_platform('CUDA')
cache.global_context_cache.platform = platform
context_cache = cache.ContextCache(capacity=None, time_to_live=None)

# Load arguments
parser = argparse.ArgumentParser(description='run repex')
parser.add_argument('out_dir', type=str, help='path to input/output dir')
parser.add_argument('sub_dir', type=str, help='sub_dir')
parser.add_argument('phase', type=str, help='phase of the simulation to use in storage filename')
parser.add_argument('replica', type=int, help='replica index')
args = parser.parse_args()

def retrieve_states(desired_state_index, iteration, replica_id):
    """
    Retrieve the thermodynamic and sampler states for a given iteration and state.

    Parameters
    ----------
    desired_state_index : int
        thermodynamic state index
    iteration : int
        iteration number
    replica_id : int
        replica id

    Returns
    -------
    thermodynamic_state : openmmtools.states.ThermodynamicState
        thermodynamic state corresponding to user specified state index and iteration
    sampler_state : openmmtools.states.SamplerState
        sampler state corresponding to user specified replica id and iteration
    """
    
    thermodynamic_state = thermodynamic_states[desired_state_index]
   
    positions = all_positions[iteration, replica_id, :, :].astype(np.float64)
    positions = unit.Quantity(positions, unit.nanometers)
    box_vectors = all_box_vectors[iteration, replica_id, :, :].astype(np.float64)
    box_vectors = unit.Quantity(box_vectors, unit.nanometers)
    sampler_state = states.SamplerState(positions=positions, box_vectors=box_vectors)

    return thermodynamic_state, sampler_state

def retrieve_potential_energy(thermodynamic_state, sampler_state, lambda_val):
    """
    Retrieve the forces for at a given lambda (for a particular configuration).

    Parameters
    ----------
    thermodynamic_state : openmmtools.states.ThermodynamicState
        thermodynamic state
    sampler_state : openmmtools.states.SamplerState
        sampler state
    lambda_val : float
        lambda value at which to retrieve the energy

    Returns
    -------
    energy : float*unit.kilojoules_per_mole
        potential energy
    """
    
    temperature = 300 * unit.kelvin
    T_max = temperature
    beta_0 = 1 / (kB * temperature)
    beta_m = 1 / (kB * T_max)
    thermodynamic_state_copy = copy.deepcopy(thermodynamic_state)
    thermodynamic_state_copy.set_alchemical_parameters(lambda_val, beta_0, beta_m)
    
    context, context_integrator = context_cache.get_context(thermodynamic_state_copy)
    sampler_state.apply_to_context(context)
    
    energy = context.getState(getEnergy=True).getPotentialEnergy()
    
    return energy

def compute_du_dlambda(u_plus_epsilon, u_minus_epsilon, epsilon=1e-3): 
    """
    Compute du/dlambda using two-sided finite difference method

    Parameters
    ----------
    u_plus_epsilon : float
        potential energy plus step size
    u_minus_epsilon : float
        potential energy minus step size
    epsilon : float, default 1e-3
        step size

    Returns
    -------
    du_dlambda : float
        du/dlambda
    """
    return(u_plus_epsilon - u_minus_epsilon) / (2 * epsilon) 

def compute_du_dlambda_onesided(u_plus_epsilon, u, epsilon=1e-3):
    """
    Compute du/dlambda using one-sided finite difference method (at lambda = 0)

    Parameters
    ----------
    u_plus_epsilon : float
        potential energy plus step size
    u : float
        potential energy
    epsilon : float, default 1e-3
        step size

    Returns
    -------
    du_dlambda : float
        du/dlambda
    """
    return (u_plus_epsilon - u) / (epsilon) 

def compute_du_dlambda_onesided_minus(u, u_minus_epsilon, epsilon=1e-3):
    """
    Compute du/dlambda using one-sided finite difference method (at lambda = 1)

    Parameters
    ----------
    u : float
        potential energy
    u_minus_epsilon : float
        potential energy minus step size
    epsilon : float, default 1e-3
        step size

    Returns
    -------
    du_dlambda : float
        du/dlambda
    """
    return (u - u_minus_epsilon) / (epsilon) 
    
def get_du_dlambda(state_index, iteration, replica_id):
    """
    Get du/dlambda given the thermodynamic state index, iteration, and replica id

    Parameters
    ----------
    desired_state_index : int
        thermodynamic state index
    iteration : int
        iteration number
    replica_id : int
        replica id

    Returns
    -------
    du_dlambda : float
        du/dlambda 
    """
    thermodynamic_state, sampler_state = retrieve_states(state_index, iteration, replica_id)
    lambda_val = lambda_schedule[state_index]
    if state_index == 0:
        u_A = retrieve_potential_energy(thermodynamic_state, sampler_state, lambda_val + epsilon)
        u_B = retrieve_potential_energy(thermodynamic_state, sampler_state, lambda_val)
        du_dlambda = compute_du_dlambda_onesided(u_A, u_B, epsilon=epsilon)
    elif state_index == n_states - 1:
        u_A = retrieve_potential_energy(thermodynamic_state, sampler_state, lambda_val)
        u_B = retrieve_potential_energy(thermodynamic_state, sampler_state, lambda_val - epsilon)
        du_dlambda = compute_du_dlambda_onesided_minus(u_A, u_B, epsilon=epsilon)
    else:
        u_A = retrieve_potential_energy(thermodynamic_state, sampler_state, lambda_val + epsilon)
        u_B = retrieve_potential_energy(thermodynamic_state, sampler_state, lambda_val - epsilon)
        du_dlambda = compute_du_dlambda(u_A, u_B, epsilon=epsilon)
    
    return du_dlambda

# Load nc files
out_dir = args.out_dir
sub_dir = args.sub_dir
phase = args.phase

reporter = MultiStateReporter(os.path.join(out_dir, f"{sub_dir}_{phase}.nc"), 'r')
thermodynamic_states = np.array(reporter.read_thermodynamic_states()[0])

nc = open_netcdf(os.path.join(out_dir, f"{sub_dir}_{phase}.nc"))
nc_checkpoint = open_netcdf(os.path.join(out_dir, f"{sub_dir}_{phase}_checkpoint.nc"))

# Get du_dlambdas
epsilon = 1e-3
n_states = nc.variables['states'].shape[1]
n_iterations = nc.variables['last_iteration'][0]
lambda_schedule = np.linspace(0.,1.,n_states)
all_positions = nc_checkpoint.variables['positions']
all_box_vectors = nc_checkpoint.variables['box_vectors']
checkpoint_interval = nc_checkpoint.CheckpointInterval
n_checkpoint_iterations = int(n_iterations / checkpoint_interval)
replica_state_indices = nc.variables['states'] # states[iteration][replica] is the thermodynamic state index (0..n_states-1) of replica 'replica' of iteration 'iteration'.

replica_index = args.replica
print(f"replica index: {replica_index}", file=sys.stderr)
du_dlambda_per_replica = []
for iteration in range(0, n_checkpoint_iterations + 1):
    print(f"iteration: {iteration}", file=sys.stderr)
    state_index = replica_state_indices[iteration*checkpoint_interval][replica_index]
    du_dlambda = get_du_dlambda(state_index, iteration, replica_index).value_in_unit_system(unit.md_unit_system)
    du_dlambda_per_replica.append(du_dlambda)
du_dlambda_per_replica = np.array(du_dlambda_per_replica)

with open(os.path.join(out_dir, f"{sub_dir}_{phase}_du_dlambda_replica{replica_index}.npy"), "wb") as f:
   np.save(f, du_dlambda_per_replica)

