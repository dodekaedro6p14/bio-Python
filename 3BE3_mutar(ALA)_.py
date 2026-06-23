from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue

# 1. Cargamos la estructura original
parser = PDBParser()
estructura = parser.get_structure('3BE3', '3BE3.pdb')

# 2. Definimos la mutación: Vamos a cambiar un residuo 
# 
# Nota: Debes verificar en tu archivo .pdb 
# el número exacto de residuo (resi)
for model in estructura:
    for chain in model:
        # Buscamos el residuo HIS en la posición 59
        if 77 in chain:
            res = chain[77]
            if res.resname == 'HIS':
                # Lo cambiamos por Alanina (ALA) - 
                # Una mutación simple
                res.resname = 'ALA'
                print("Mutación realizada: HIS77 -> ALA77")

# 3. Guardamos la nueva estructura para abrirla en PyMOL
io = PDBIO()
io.set_structure(estructura)
io.save("3BE3_mutado.pdb")