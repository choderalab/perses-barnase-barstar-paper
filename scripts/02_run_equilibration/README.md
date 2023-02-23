# Run equilibration
Run gentle equilibration based on the protocol outlined in Table S2 [here](https://academic.oup.com/nar/article/44/1/63/2499624?login=true#supplementary-data).
Outputs equilibrated PDB file and the corresponding box vectors as a .npy file.

These files are then used as input for generating the input files for free energy calculations (Perses `RESTCapableHybridTopologyFactory` objects).

## Python scripts
- `run_equilibration.py` - runs equilibration

## Bash scripts
Bash scripts for running equilibration for each experiment are located in `perses-barnase-barstar-paper/data/`. 
