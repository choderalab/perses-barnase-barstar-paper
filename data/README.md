# Data

Note "htf" refers to the Perses `RESTCapableHybridTopologyFactory` object, which contains the hybrid topology, positions, and system for a given an alchemical transformation.

## File organization
Bash scripts for running experiments and data for generating figures are organized according to this structure: `{batch_number}/{experiment_number}/{replicate_number}/`.
Batches, experiments, and replicates are described in detail below.

Tables of $\Delta\Delta G$ values for terminally-blocked amino acid and barnase:barstar mutations are available in:
- `table_10ns_arex.csv` -- barnase:barstar mutations (10 ns/replica AREX for complex phase, 10 ns/replica AREX for apo phase)
- `table_10ns_arest.csv` -- barnase:barstar mutations (10 ns/replica AREST for complex phase, 10 ns/replica AREX for apo phase)
- `table_50ns_arest.csv` -- barnase:barstar mutations (50 ns/replica AREST for complex phase, 10 ns/replica AREX for apo phase)
- `table_50ns_arex.csv` -- barnase:barstar mutations (50 ns/replica AREX for complex phase, 10 ns/replica AREX for apo phase)
- `table_terminally_blocked.csv` -- terminally-blocked amino acid mutations (5 ns/replica AREX for both ACE-X-NME and ALA-X-ALA phases)

Tables of the names of the degree of freedoms corresponding to PCCs are available in:
- `table_10ns_arex_pccs.csv` -- barnase:barstar mutations (10 ns/replica AREX for complex phase)
- `table_50ns_arest_pccs.csv` -- barnase:barstar mutations (50 ns/replica AREST for complex phase)

Details on the `Degree of freedom` column: For backbone and sidechain torsions, the degree of freedom is formatted according to "{chain name}-{residue id} {dihedral angle name} {'old' or 'new', aka WT and mutant residues, respectively}", where the dihedral angle name can be phi, psi, chi1, chi2, chi3, or chi4. For intra and inter interface contacts, the degree of freedom is formatted according to "{chain name of residue A}-{residue A id} | {chain name of residue B}-{residue B id} {'old' or 'new', aka WT and mutant topologies, respectively}". For neighboring waters, the degree of freedom is formatted according to "nearby waters {'old' or 'new', aka WT and mutant residues, respectively}", where "nearby waters" refers to waters within 5 angstroms of the WT or mutant residue. 

Table of $\Delta G$ values for computing the multi-protonation state $\Delta\Delta G$ s for barnase:barstar mutations D35A and K27A is available at: `D35A_K27A.xlsx`.

## Batches
- `43` - Contains scripts for solvating/parametrizing and running equilibration
- `45` - Contains for scripts for generating htfs, running/resuming AREX, analyzing free energies, and running $\partial U$ / $\partial \lambda$ analysis for barnase:barstar mutations A42T, R87A, and W44F. Also contains pickle files for free energy data and $\partial U$ / $\partial \lambda$ correlation data.
- `46` - Contains for scripts for generating htfs, running AREX, and analyzing free energies for all terminally-blocked amino acid mutations. Also contains pickle files for free energy data.
- `47` - Contains for scripts for generating htfs, running/resuming AREX, analyzing free energies, and running $\partial U$ / $\partial \lambda$ analysis for all barnase:barstar mutations besides A42T, R87A, and W44F. Also contains pickle files for free energy data and $\partial U$ / $\partial \lambda$ correlation data.
- `48` - Contains for scripts for running AREX (with restraints) and analyzing free energies for barnase:barstar mutations A42T and R87A. Also contains pickle files for free energy data. Htfs were copied from `45/9/replicate_1` and `45/10/replicate_1` for `48/0/replicate_0` and `48/1/replicate_1`, respectively.
- `49` - Contains for scripts for generating modified htfs, running AREST, and analyzing free energies for barnase:barstar mutations A42T and R87A across 9 REST parameter combinations for each mutation. Also contains pickle files for free energy data. Htfs were copied from `45/9/replicate_1` for `49/0-2/replicate_0-2` and `45/10/replicate_1` for `49/9-11/replicate_0-2`. Htfs were also copied from replicate 0 to replicates 1 and 2.
- `52` - Contains for scripts for generating htfs, running/resuming AREST, and analyzing free energies for all barnase:barstar mutations. Also contains pickle files for free energy data. Htfs were copied from `49/4/replicate_0` and `49/13/replicate_0` for `52/12/replicate_0` and `52/17/replicate_0`, respectively.

## Experiments
- `43`
    - `0` - ACE-ALA-NME in solvent
	- `1` - ACE-ARG-NME in solvent
	- `2` - ACE-TRP-NME in solvent
	- `3` - WT barstar and barstar:barnase
	- `4` - WT barnase and barnase:barstar
	- `5` - A42T barstar and A42T barstar:barnase
	- `6` - ACE-TYR-NME in solvent
	- `7` - ACE-THR-NME in solvent
	- `8` - ACE-GLU-NME in solvent
	- `9` - ACE-ASP-NME in solvent
	- `10` - ACE-LYS-NME in solvent
	- `11` - ACE-HIS-NME in solvent
	- `12` - ACE-PHE-NME in solvent
	- `13` - ACE-GLN-NME in solvent
	- `14` - ALA-ALA-ALA in solvent
	- `15` - ALA-ARG-ALA in solvent
	- `16` - ALA-TRP-ALA in solvent
	- `17` - ALA-TYR-ALA in solvent
	- `18` - ALA-THR-ALA in solvent
	- `19` - ALA-GLU-ALA in solvent
	- `20` - ALA-ASP-ALA in solvent
	- `21` - ALA-LYS-ALA in solvent
	- `22` - ALA-HIS-ALA in solvent
	- `23` - ALA-PHE-ALA in solvent
	- `24` - ALA-GLN-ALA in solvent
	- `25` - F29 barstar and F29 barstar:barnase
	- `26` - A29 barstar and A29 barstar:barnase
	- `27` - F38 barstar and F38 barstar:barnase
	- `28` - F44 barstar and F44 barstar:barnase
	- `29` - A76 barstar and A76 barstar:barnase
	- `30` - A80 barstar and A80 barstar:barnase
	- `31` - A35 barstar and A35 barstar:barnase
	- `32` - A39 barstar and A39 barstar:barnase
	- `33` - A27 barnase and A27 barnase:barstar
	- `34` - A59 barnase and A59 barnase:barstar
	- `35` - Q83 barnase and Q83 barnase:barstar
	- `36` - A87 barnase and A87 barnase:barstar
	- `37` - A102 barnase and A102 barnase:barstar
	- `38` - ASH capped peptide  in solvent
	- `40` - ASH35 barstar and ASH35 barstar:barnase
	- `43` - LYN capped peptide in solvent
	- `44` - LYN27 barnase and LYN27 barnase:barstar

- `45`
    - `9` - A42T barstar, A42T barstar:barnase
    - `10` - R87A barnase, R87A barnase:barstar
    - `11` - W44F barstar, W44F barstar:barnase

- `46`
    - (ACE-X-NME)
        - `0` - Y->F
        - `1` - Y->A
        - `2` - W->F
        - `3` - T->A
        - `4` - E->A
        - `5` - D->A
        - `6` - K->A
        - `7` - R->A
        - `8` - R->Q
        - `9` - H->A
        - `10` - F->Y
        - `11` - A->Y
        - `12` - F->W
        - `13` - A->T
        - `14` - A->E
        - `15` - A->D
        - `16` - A->K
        - `17` - A->R
        - `18` - Q->R
        - `19` - A->H 
        - `40` - ASH->A
        - `44` - LYN->A
    - (ALA-X-ALA)
        - `20` - Y->F
        - `21` - Y->A
        - `22` - W->F
        - `23` - T->A
        - `24` - E->A
        - `25` - D->A
        - `26` - K->A
        - `27` - R->A
        - `28` - R->Q
        - `29` - H->A
        - `30` - F->Y
        - `31` - A->Y
        - `32` - F->W
        - `33` - A->T
        - `34` - A->E
        - `35` - A->D
        - `36` - A->K
        - `37` - A->R
        - `38` - Q->R
        - `39` - A->H 
- `47`
	- `0` - Y29F
	- `1` - Y29A
	- `2` - W38F
	- `3` - T42A
	- `4` - W44F -- note this directory is empty, see `45/11`
	- `5` - E76A
	- `6` - E80A
	- `7` - D35A
	- `8` - D39A
	- `9` - K27A
	- `10` - R59A
	- `11` - R83Q
	- `12` - R87A -- note this directory is empty, see `45/10`
	- `13` - H102A
	- `14` - F29Y
	- `15` - A29Y
	- `16` - F38W
	- `17` - A42T -- note this directory is empty, see `45/9`
	- `18` - F44W
	- `19` - A76E
	- `20` - A80E
	- `21` - A35D
	- `22` - A39D
	- `23` - A27K
	- `24` - A59R
	- `25` - Q83R
	- `26` - A87R
	- `27` - A102H
	- `28` - ASH35A
	- `32` - LYN27A
- `48`
	- `0` - A42T 
	- `1` - R87A
- `49`
	- `0` - A42T, radius 0.3 nm, T_max 400 K
	- `1` - A42T, radius 0.3 nm, T_max 600 K
	- `2` - A42T, radius 0.3 nm, T_max 1200 K
	- `3` - A42T, radius 0.5 nm, T_max 400 K
	- `4` - A42T, radius 0.5 nm, T_max 600 K
	- `5` - A42T, radius 0.5 nm, T_max 1200 K
	- `6` - A42T, radius 0.7 nm, T_max 400 K
	- `7` - A42T, radius 0.7 nm, T_max 600 K
	- `8` - A42T, radius 0.7 nm, T_max 1200 K
	- `9` - R87A, radius 0.3 nm, T_max 400 K
	- `10` - R87A, radius 0.3 nm, T_max 600 K
	- `11` - R87A, radius 0.3 nm, T_max 1200 K
	- `12` - R87A, radius 0.5 nm, T_max 400 K
	- `13` - R87A, radius 0.5 nm, T_max 600 K
	- `14` - R87A, radius 0.5 nm, T_max 1200 K
	- `15` - R87A, radius 0.7 nm, T_max 400 K
	- `16` - R87A, radius 0.7 nm, T_max 600 K
	- `17` - R87A, radius 0.7 nm, T_max 1200 K
- `52`
	- `0` - Y29F
	- `1` - Y29A
	- `2` - W38F
	- `3` - T42A
	- `4` - W44F
	- `5` - E76A
	- `6` - E80A
	- `7` - D35A
	- `8` - D39A
	- `9` - K27A
	- `10` - R59A
	- `11` - R83Q
	- `12` - R87A
	- `13` - H102A
	- `14` - F29Y
	- `15` - A29Y
	- `16` - F38W
	- `17` - A42T
	- `18` - F44W
	- `19` - A76E
	- `20` - A80E
	- `21` - A35D
	- `22` - A39D
	- `23` - A27K
	- `24` - A59R
	- `25` - Q83R
	- `26` - A87R
	- `27` - A102H

## Replicates
Below, we list the names of the replicate directories used for each batch.
- `45` - replicate 1 for `9` and `10`, replicate 0 for `11`
- `46` - replicate 0
- `47` - replicate 1
- `48` - replicate 0 for `0` and replicate 1 for `1`
- `49` - replicates 0, 1, 2
- `52` - replicate 0 for all except replicate 1 for `19` and `21` 
