import os
import shutil
import tempfile
import logging
from openmm import app
from openmoltools.utils import getoutput

# Note the below function is modified from openmoltools.amber.build_peptide_tleap
def build_peptide_tleap(amino_acids):
    """
    Use tleap to generate a peptide topology
    and positions.
    Parameters
    ----------
    amino_acids : list of str
        List of amino acids and caps in three-letter names
        e.g. ['ACE', 'ALA', 'NME']
    Returns
    -------
    topology : simtk.openmm.app.Topology object
        Topology of the amino acid
    positions : [n, 3] array
        positions of atoms
    """
    aa_str = " ".join(amino_acids)
    filename = "".join(amino_acids)
    tleapstr = """
    source oldff/leaprc.ff99SBildn
    system = sequence {{ {amino_acid} }}
    impose system {{1 2 3}} {{{{"N" "CA" "C" "N" -40.0}} {{"C" "N" "CA" "C" -60.0}}}} # Impose alpha helical conformation
    saveamberparm system {filename}.prmtop {filename}.inpcrd
    """.format(amino_acid=aa_str, filename=filename)
    cwd = os.getcwd()
    temp_dir = tempfile.mkdtemp()
    os.chdir(temp_dir)
    tleap_file = open('tleap_commands', 'w')
    tleap_file.writelines(tleapstr)
    tleap_file.close()
    tleap_cmd_str = "tleap -f %s " % tleap_file.name

    #call tleap, log output to logger
    output = getoutput(tleap_cmd_str)
    logging.debug(output)

    prmtop = app.AmberPrmtopFile("{filename}.prmtop".format(filename=filename))
    inpcrd = app.AmberInpcrdFile("{filename}.inpcrd".format(filename=filename))
    topology = prmtop.topology
    positions = inpcrd.positions

    os.chdir(cwd)
    shutil.rmtree(temp_dir)
    return topology, positions


amino_acids = ['ALA', 'ARG', 'TRP', 'TYR', 'THR', 'GLU', 'ASP', 'LYS', 'HIS', 'PHE', 'GLN']
for amino_acid in amino_acids:
    for i, sequence in enumerate([["ACE", amino_acid, "NME"], ["NALA", amino_acid, "CALA"]]):
        print(f"generating {sequence}")       
        sequence_type = '' if i == 0 else '_ala'
        topology, positions = build_peptide_tleap(sequence)
        app.PDBFile.writeFile(topology, positions, open(f"{amino_acid}_capped{sequence_type}.pdb", "w"), keepIds=True)
