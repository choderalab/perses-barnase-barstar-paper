# Generate htfs
Generate the input files for free energy calculations (Perses `RESTCapableHybridTopologyFactory` objects aka rest-capable htfs) and save as pickle files. 

These files are then used as input for running repex simulations.

## Python scripts
- generate_htfs.py - generate rest-capable htfs
- generate_modified_htfs.py - take an existing rest-capable htf's topology_proposal and generate a rest-capable htf with the specified rest radius and/or w_lifting modified
- generate_vanilla_htfs.py - take an existing rest-capable htf's topology_proposal and generate a vanilla htf (Perses `HybridTopologyFactory`)

## Bash scripts
Bash scripts for generating the htfs for each experiment are located in `perses-barnase-barstar-paper/data/`.
