# Input files

## Files
- `1brs_barnase_{mutant}.pdb` - Mutant barnase, where `mutant` is the one letter amino acid code followed by the residue position (e.g. A102 for alanine at position 102)
- `1brs_barnase_renumbered.pdb` - WT barnase
- `1brs_barstar_{mutant}.pdb` - Mutant barstar, where `mutant` is the one letter amino acid code followed by the residue position (e.g. A42 for alanine at position 42)
- `1brs_barstar_renumbered.pdb` - WT barstar
- `{amino_acid_name}_capped.pdb` - terminally blocked amino acids with ACE and NME caps where `amino_acid_name` is the three letter amino acid code for the amino acid of interest (e.g. TRP for tryptophan). For TRP, the amino acid sequence is: ACE-TRP-NME
- `{amino_acid_name}_capped_ala.pdb` - terminally blocked amino acids with ALA caps where `amino_acid_name` is the three letter amino acid code for the amino acid of interest (e.g. TRP for tryptophan). For TRP, the aminon acid sequence is: ALA-TRP-ALA
- `generate_nonstandard_protonation_states.py` - use openmm to generate ASH and LYN PDBs
- `generate_peptide_pdbs.py` - use tleap to generate the ACE-X-NME and ALA-X-ALA pdbs
- `renumber.py` - renumber the residues starting at the NME cap in the maestro prepped WT barnase and barstar PDBs

## Structure preparation workflow
### Barnase:barstar

Preparation of WT structures (in Schrodinger Maestro 2021-2)
1. Load PDB ID 1BRS
2. Keep chains A and D by manually selecting other chains and deleting them
3. Run protein prep wizard:
   - Add caps
   - Protonate
     - HIS18 (in barnase) was protonated as HID, HIS102 (in barnase) was protonated as HIE, and H17 (in barstar) was protonated as HID.
   - Add in missing loop (64-65 of chain D)
   - Optimize hydrogen bond network and orientations of ASN, GLN, and HIS residues at pH 8.0
4. Split structure into barnase and barstar
	- Saved as 1brs_barnase.pdb (chain A) and 1brs_barstar.pdb (chain D)
	- Ran renumber.py to renumber NME residues to generate 1brs_barnase_renumbered.pdb and 1brs_barstar_renumberd.pdb

Preparation of mutant structures (in Schrodinger Maestro 2021-3)
1. Load 1brs_barnase_renumbered.pdb (or 1brs_barstar_renumbered.pdb)
2. Right click on residue to mutate and click "Mutate Residue", then choose residue to mutate to
3. Make sure that the sidechain lines up well with WT sidechain and if not, right click and choose "Select Rotamer", then select rotamer that most closely lines up with WT sidechain
4. Compare the mutant rotamer to existing crystal structures, if available  
5. Save PDB

### Peptide
```
python generate_peptide_pdbs.py
```
