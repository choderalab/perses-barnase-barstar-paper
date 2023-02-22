# Input files

## Files
- 1brs_barnase_{mutant}.pdb - Mutant barnase, where `mutant` is the one letter amino acid code followed by the residue position (e.g. A102 for alanine at position 102)
- 1brs_barnase_renumbered.pdb - WT barnase
- 1brs_barstar_{mutant}.pdb - Mutant barstar, where `mutant` is the one letter amino acid code followed by the residue position (e.g. A42 for alanine at position 42)
- 1brs_barstar_renumbered.pdb - WT barstar
- {amino_acid_name}_capped.pdb - terminally blocked amino acids with ACE and NME caps where `amino_acid_name` is the three letter amino acid code for the amino acid of interest (e.g. TRP for tryptophan). For TRP, the amino acid sequence is: ACE-TRP-NME
- {amino_acid_name}_capped_ala.pdb - terminally blocked amino acids with ALA caps where `amino_acid_name` is the three letter amino acid code for the amino acid of interest (e.g. TRP for tryptophan). For TRP, the aminon acid sequence is: ALA-TRP-ALA
- generate_peptide_pdbs.py - use tleap to generate the ACE-X-NME and ALA-X-ALA pdbs
- renumber.py - renumber the residues starting at the NME cap in the maestro prepped WT barnase and barstar PDBs

## Structure preparation workflow
### Barnase:barstar
Note that I think I fed perses_protein_mutations/input/1brs_{barstar/barnase}.pdb to renumber.py to get the current files.
### Peptide

