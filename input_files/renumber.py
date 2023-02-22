import argparse
from openmm import app

# Read args
parser = argparse.ArgumentParser(description='renumber PDB starting from NME and onwards because maestro add caps as inserted residues')
parser.add_argument('input_file', type=str, help='path to input file')
args = parser.parse_args()

pdb = app.PDBFile(args.input_file)
for residue in pdb.topology.residues():
    if residue.name == 'NME' or residue.name == 'HOH':
        print(residue.name)
        residue.id = str(int(residue.id) + 1)
        residue.insertionCode = ' '
app.PDBFile.writeFile(pdb.topology, pdb.positions, open(args.input_file[:-4] + "_renumbered.pdb", "w"), keepIds=True)

