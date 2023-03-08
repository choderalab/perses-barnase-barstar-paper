# perses-barnase-barstar-paper
Contains scripts, input files, and data for:
- Running and analyzing relative free energy calculations for terminally-blocked amino acids and barnase:barstar using [Perses](https://github.com/choderalab/perses)
- Identifying and characterizing slow degrees of freedom

Preprint: https://doi.org/10.1101/2023.03.07.530278

## Contributors
- Ivy Zhang
- Dominic Rufa

## License
* This software is licensed under the [MIT license](https://opensource.org/licenses/MIT) - a copy of this license is provided as `SOFTWARE_LICENSE`
* The data in this repository is made available under the Creative Commons [CC0 (“No Rights Reserved”) License](https://creativecommons.org/share-your-work/public-domain/cc0/) - a copy of this license is provided as `DATA_LICENSE`

## Manifest
* `data` - Contains the data necessary to reproduce the figures and bash scripts used to run the jobs
* `envs` - Contains the dumps of the conda environments used in this study
* `input_files` - Contains all input PDB files and relevant scripts for generating/modifying the PDBs
* `scripts` - Contains all Python scripts for solvating, parametrizing, running equilibration, generating free energy calculation input files, running alchemical replica exchange (AREX), and running alchemical replica exchange with solute tempering (AREST).

## Environment
Core dependencies include Perses 0.10.1, OpenMMTools 0.21.5, MDTraj 1.9.7, and pymbar 3.1.1.

[OpenMM 8.0.0beta](https://anaconda.org/conda-forge/openmm/files?version=8.0.0beta) (build 0), a development version of OpenMM 7, was used to generate the input files for alchemical replica exchange (AREX) and alchemical replica exchange with solute tempering (AREST), run equilibration, and run AREX for the terminally-blocked amino acids. The conda environment including this version of OpenMM is called `perses-paper3`.

[OpenMM 7.7.0.dev2](https://anaconda.org/conda-forge/openmm/files?version=7.7.0dev2), a development version of OpenMM 7 which was built after OpenMM 8.0.0beta and contains a performance enhancement for AREX and AREST, was used for all other AREX and AREST simulations. The conda environment including this version of OpenMM is called `perses-paper5`. This environment was also used to conduct all analysis.

$\Delta\Delta G$ comparison plots were generated with cinnabar 0.3.0. 
All other plots were generated using Matplotlib 3.5.2.

The environment files for `perses-paper3` and `perses-paper5` are located in `envs`.
They can be used to create a new environment by:
```
mamba create -n new-env
mamba install --name new-env --file perses-paper5-explicit.txt
```
