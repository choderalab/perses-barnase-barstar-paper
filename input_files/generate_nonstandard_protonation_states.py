from openmm import app, unit

## ASH ##

files = ["1brs_barstar_renumbered.pdb", "ASP_capped.pdb"]
resids = ["35", "2"]

for f, resid in zip(files, resids):
    # Load in WT 
    pdb = app.PDBFile(f)

# Create list of variants to use for each residue. None means the variant will be selected automatically
    variants = []
    for res in pdb.topology.residues():
        if res.id == resid and res.name == 'ASP':
            variants.append('ASH')
        else:
            variants.append(None)

    # Add hydrogens with the above variants list to protonate D35 appropriately (s.t. its ASH)
    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(pH=8.0, variants=variants)

    # Save
    outfile = "1brs_barstar_ASH35.pdb" if "barstar" in f else ("ASH_capped.pdb" if "ala" not in f else "ASH_capped_ala.pdb") 
    app.PDBFile.writeFile(modeller.topology, modeller.positions, open(outfile, "w"), keepIds=True)

## LYN ##

files = ["1brs_barnase_renumbered.pdb", "LYS_capped.pdb"]
resids = ["27", "2"]

for f, resid in zip(files, resids):
    # Load in WT
    pdb = app.PDBFile(f)

    # Create list of variants to use for each residue. None means the variant will be selected automatically
    variants = []
    for res in pdb.topology.residues():
        if res.id == resid and res.name == 'LYS':
            variants.append('LYN')
        else:
            variants.append(None)

    # Add hydrogens with the above variants list to deprotonate K27 appropriately (s.t. its LYN)
    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(pH=8.0, variants=variants)

    # Save
    outfile = "1brs_barnase_LYN27.pdb" if "barnase" in f else ("LYN_capped.pdb" if "ala" not in f else "LYN_capped_ala.pdb")
    app.PDBFile.writeFile(modeller.topology, modeller.positions, open(outfile, "w"), keepIds=True)
