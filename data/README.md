# Data

## File organization
`{batch_number}/{experiment_number}/{replicate_number}/`

## Batches
- `43` - Contains scripts for solvating/parametrizing and running equilibration
- `45` - Contains for scripts for generating htfs, running/resuming AREX, analyzing free energies, and running $\partial U$ / $\partial \lambda$ analysis for barnase:barstar mutations A42T, R87A, and W44F. Also contains pickle files for free energy data and $\partial U$ / $\partial \lambda$ correlation data.
- `46` - Contains for scripts for generating htfs, running AREX, and analyzing free energies for all terminally-blocked amino acid mutations. Also contains pickle files for free energy data.
- `47` - Contains for scripts for generating htfs, running/resuming AREX, analyzing free energies, and running $\partial U$ / $\partial \lambda$ analysis for all barnase:barstar mutations besides A42T, R87A, and W44F. Also contains pickle files for free energy data and $\partial U$ / $\partial \lambda$ correlation data.
- `48` - Contains for scripts for running AREX (with restraints) and analyzing free energies for barnase:barstar mutations A42T and R87A. Also contains pickle files for free energy data.
- `49` - Contains for scripts for running AREST and analyzing free energies for barnase:barstar mutations A42T and R87A across 9 REST parameter combinations for each mutation. Also contains pickle files for free energy data.
- `52` - Contains for scripts for generating htfs, running/resuming AREST, and analyzing free energies for all barnase:barstar mutations. Also contains pickle files for free energy data.

## Experiments
- `43`
    - 0 - ACE-ALA-NME in solvent
    - 1 - ACE-ARG-NME in solvent
    - 2 - ACE-TRP-NME in solvent
    - 3 - WT barstar and barstar:barnase
    - 4 - WT barnase and barnase:barstar
    - 5 - A42T barstar and A42T barstar:barnase
    - 6 - ACE-TYR-NME in solvent
    - 7 - ACE-THR-NME in solvent
    - 8 - ACE-GLU-NME in solvent
    - 9 - ACE-ASP-NME in solvent
    - 10 - ACE-LYS-NME in solvent
    - 11 - ACE-HIS-NME in solvent
    - 12 - ACE-PHE-NME in solvent
    - 13 - ACE-GLN-NME in solvent
    - 14 - ALA-ALA-ALA in solvent
    - 15 - ALA-ARG-ALA in solvent
    - 16 - ALA-TRP-ALA in solvent
    - 17 - ALA-TYR-ALA in solvent
    - 18 - ALA-THR-ALA in solvent
    - 19 - ALA-GLU-ALA in solvent
    - 20 - ALA-ASP-ALA in solvent
    - 21 - ALA-LYS-ALA in solvent
    - 22 - ALA-HIS-ALA in solvent
    - 23 - ALA-PHE-ALA in solvent
    - 24 - ALA-GLN-ALA in solvent
    - 25 - F29 barstar and F29 barstar:barnase
    - 26 - A29 barstar and A29 barstar:barnase
    - 27 - F38 barstar and F38 barstar:barnase
    - 28 - F44 barstar and F44 barstar:barnase
    - 29 - A76 barstar and A76 barstar:barnase
    - 30 - A80 barstar and A80 barstar:barnase
    - 31 - A35 barstar and A35 barstar:barnase
    - 32 - A39 barstar and A39 barstar:barnase
    - 33 - A27 barnase and A27 barnase:barstar
    - 34 - A59 barnase and A59 barnase:barstar
    - 35 - Q83 barnase and Q83 barnase:barstar
    - 36 - A87 barnase and A87 barnase:barstar
    - 37 - A102 barnase and A102 barnase:barstar
    - 38 - ASH capped peptide  in solvent
    - 40 - ASH35 barstar and ASH35 barstar:barnase
    - 43 - LYN capped peptide in solvent
    - 44 - LYN27 barnase and LYN27 barnase:barstar

## Replicates
Below, we list the names of the replicate directories used for each batch.
- `45` - replicate 1 for `9` and `10`, replicate 0 for `11`
- `46` - replicate 0
- `47` - replicate 1
- `48` - replicate 0
- `49` - replicates 0, 1, 2
- `52` - replicate 0 for all except replicate 1 for `19` and `21` 
