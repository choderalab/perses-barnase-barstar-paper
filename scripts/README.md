# Scripts
The scripts are organized based on each step in the workflow:
1. `01_generate_solvated_inputs` - Generate PDB and OpenMM system.xml files which contain OpenMM topology, positions, and system for a WT system of interest
2. `02_run_equilibration` - Run gentle equilibration
3. `03_generate_htfs` - Generate Perses RESTCapableHybridTopologyFactory objects which contain the hybrid topology, positions, and system
4. `04_run_repex` - Run replica exchange simulations (and resume, if necessary)
5. `05_analyze` -  Analyze the simulations (compute free energies, compute $\partial U$ / $\partial \lambda$, monitor degrees of freedom, check replica mixing, check slopes of free energy time series, etc.)

