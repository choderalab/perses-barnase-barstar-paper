# perses-barnase-barstar-paper
Contains scripts, input files, and data for:
- Running and analyzing relative free energy calculations for terminally-blocked amino acids and barnase:barstar using [Perses](https://github.com/choderalab/perses)
- Identifying and characterizing slow degrees of freedom

Publication: Coming soon!

## Contributors
- Ivy Zhang
- Dominic Rufa

## License
* This software is licensed under the [MIT license](https://opensource.org/licenses/MIT) - a copy of this license is provided as `SOFTWARE_LICENSE`
* The data in this repository is made available under the Creative Commons [CC0 (“No Rights Reserved”) License](https://creativecommons.org/share-your-work/public-domain/cc0/) - a copy of this license is provided as `DATA_LICENSE`

## Manifest
* `data` - Contains the data necessary to reproduce the figures
* `input_files` - Contains all input PDB files and relevant scripts for generating/modifying the PDBs
* `scripts` - Contains all scripts for solvating, parametrizing, running equilibration, generating free energy calculation input files, running alchemical replica exchange (AREX), and running alchemical replica exchange with solute tempering (AREST).

## Environment
Core dependencies include Perses 0.10.1, OpenMMTools 0.21.5, MDTraj 1.9.7, and pymbar 3.1.1.

[OpenMM 8.0.0beta](https://anaconda.org/conda-forge/openmm/files?version=8.0.0beta) (build 0), a development version of OpenMM 7, was used to generate the input files for alchemical replica exchange (AREX) and alchemical replica exchange with solute tempering (AREST), run equilibration, and run AREX for the terminally-blocked amino acids. 

[OpenMM 7.7.0.dev2](https://anaconda.org/conda-forge/openmm/files?version=7.7.0dev2), a development version of OpenMM 7 which was built after OpenMM 8.0.0beta and contains a performance enhancement for AREX and AREST, was used for all other AREX and AREST simulations.

$\Delta\Delta G$ comparison plots were generated with cinnabar 0.3.0. 
All other plots were generated using Matplotlib 3.5.2.

Commands used to generate the environment:
```bash

```

The environment yaml is located at XX.
