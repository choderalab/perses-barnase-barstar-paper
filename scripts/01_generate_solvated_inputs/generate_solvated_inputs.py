import os
import argparse
import mdtraj as md
import numpy as np
from pathlib import Path
import openmm
from openmm import unit, app
from openmmforcefields.generators import SystemGenerator

# Read args
parser = argparse.ArgumentParser(description='generate solvated topology, positions, system')
parser.add_argument('protein_filename', type=str, help='path to protein file')
parser.add_argument('outdir', type=str, help='output directory')
parser.add_argument('--ligand_input', type=str, help='path to ligand_input file')
parser.add_argument('--padding', type=float, help='solvent padding')
args = parser.parse_args()

# Load PDBs
protein_pdb = app.PDBFile(args.protein_filename)
protein_positions, protein_topology, protein_md_topology = protein_pdb.positions, protein_pdb.topology, md.Topology.from_openmm(protein_pdb.topology)
protein_topology = protein_md_topology.to_openmm()
protein_n_atoms = protein_md_topology.n_atoms
if Path(args.protein_filename).stem[-3:] == 'ala': # If the pdb has ALA caps, change cap names
    residues = list(protein_topology.residues())
    residues[0].name = 'NALA'
    residues[-1].name = 'CALA'

if args.ligand_input is not None:
    ligand_pdb = app.PDBFile(args.ligand_input)
    ligand_positions, ligand_topology, ligand_md_topology = ligand_pdb.positions, ligand_pdb.topology, md.Topology.from_openmm(ligand_pdb.topology)
    ligand_n_atoms = ligand_md_topology.n_atoms

    # Combine topologies and positions
    complex_md_topology = protein_md_topology.join(ligand_md_topology)
    complex_topology = complex_md_topology.to_openmm()
    complex_positions = unit.Quantity(np.zeros([protein_n_atoms + ligand_n_atoms, 3]), unit=unit.nanometers)
    complex_positions[:protein_n_atoms, :] = protein_positions
    complex_positions[protein_n_atoms:, :] = ligand_positions

    # Convert positions back to openmm vec3 objects
    complex_positions_vec3 = []
    for position in complex_positions:
        complex_positions_vec3.append(openmm.Vec3(*position.value_in_unit_system(unit.md_unit_system)))
    complex_positions = unit.Quantity(value=complex_positions_vec3, unit=unit.nanometer)

# Create a system generator
forcefield_files = ['amber14/protein.ff14SB.xml', 'amber14/tip3p.xml']
barostat = openmm.MonteCarloBarostat(1.0 * unit.atmosphere, 300 * unit.kelvin, 50)
forcefield_kwargs = {'removeCMMotion': False, 'constraints' : app.HBonds, 'rigidWater': True, 'hydrogenMass' : 3 * unit.amus}
periodic_forcefield_kwargs = {'nonbondedMethod': app.PME, 'ewaldErrorTolerance': 0.00025}
nonperiodic_forcefield_kwargs = None
system_generator = SystemGenerator(forcefields=forcefield_files,
                                                barostat=barostat,
                                                forcefield_kwargs=forcefield_kwargs,
                                                periodic_forcefield_kwargs=periodic_forcefield_kwargs,
                                                nonperiodic_forcefield_kwargs=nonperiodic_forcefield_kwargs)

# Set solvation parameters
water_model = 'tip3p'
padding = 1.2 * unit.nanometers if args.padding is None else args.padding * unit.nanometers
print(f"Using padding: {padding}")
box_shape = 'cube'
ionic_strength = 0.05 * unit.molar

# Solvate and generate system for apo and then complex (if ligand_input is specified)
topologies_all = [protein_topology]
positions_all = [protein_positions]
phases_all = ['apo']
if args.ligand_input is not None:
    topologies_all.append(complex_topology)
    positions_all.append(complex_positions)
    phases_all.append('complex')

for topology, positions, phase in zip(topologies_all, positions_all, phases_all):
    modeller = app.Modeller(topology, positions)
    modeller.addSolvent(system_generator.forcefield, model=water_model, padding=padding, boxShape=box_shape, ionicStrength=ionic_strength)
    solvated_topology = modeller.getTopology()
    solvated_positions = modeller.getPositions()

    # Canonicalize the solvated positions: turn tuples into np.array
    solvated_positions = unit.quantity.Quantity(value=np.array([list(atom_pos) for atom_pos in solvated_positions.value_in_unit_system(unit.md_unit_system)]), unit=unit.nanometers)

    # Generate system
    solvated_system = system_generator.create_system(solvated_topology)

    # Fix naked charges
    atom_map = {atom.index: atom.residue.name for atom in solvated_topology.atoms()}
    force_dict = {force.__class__.__name__: force for force in solvated_system.getForces()}
    if 'NonbondedForce' in [k for k in force_dict.keys()]:
        nb_force = force_dict['NonbondedForce']
        for idx in range(nb_force.getNumParticles()):
            if atom_map[idx] in ['HOH', 'WAT']: # Do not add naked charge fix to water hydrogens
                continue
            charge, sigma, epsilon = nb_force.getParticleParameters(idx)
            if sigma == 0*unit.nanometer:
                new_sigma = 0.06*unit.nanometer
                nb_force.setParticleParameters(idx, charge, new_sigma, epsilon)
                print(f"Changed particle {idx}'s sigma from {sigma} to {new_sigma}")
            if epsilon == 0*unit.kilojoule_per_mole:
                new_epsilon = 0.0001*unit.kilojoule_per_mole
                nb_force.setParticleParameters(idx, charge, sigma, new_epsilon)
                print(f"Changed particle {idx}'s epsilon from {epsilon} to {new_epsilon}")
                if sigma == 1.0 * unit.nanometer: # in protein.ff14SB, hydroxyl hydrogens have sigma=1 and epsilon=0
                    new_sigma = 0.1*unit.nanometer
                    nb_force.setParticleParameters(idx, charge, new_sigma, epsilon)
                    print(f"Changed particle {idx}'s sigma from {sigma} to {new_sigma}")

    if Path(args.protein_filename).stem[-3:] == 'ala': # If the pdb has ALA caps, change names back
        for res in solvated_topology.residues():
            if res.name in ['NALA', 'CALA']:
                res.name = 'ALA'

    with open(os.path.join(args.outdir, f"{phase}_solvated.pdb"), "w") as f:
        app.PDBFile.writeFile(solvated_topology, solvated_positions, f, keepIds=True)

    with open(os.path.join(args.outdir, f"{phase}_solvated.xml"), "w") as f:
        xml = openmm.XmlSerializer.serialize(solvated_system)
        f.write(xml)
