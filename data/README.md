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
TODO: describe what each experiment corresponds to

## Replicates
Below, we list the names of the replicate directories used for each batch.
- `45` - replicate 1 for `9` and `10`, replicate 0 for `11`
- `46` - replicate 0
- `47` - replicate 1
- `48` - replicate 0
- `49` - replicates 0, 1, 2
- `52` - replicate 0 for all except replicate 1 for `19` and `21` 
