from Bio.PDB import PDBParser, PDBIO

parser = PDBParser()
estructura = parser.get_structure('3BE3', '3BE3.pdb')

# Acceder directamente a la cadena A, residuo 77
residuo_a_mutar = estructura[0]['A'][77]
residuo_a_mutar.resname = 'TRP' # Cambiamos HIS por TRP 
#(Triptófano, el aminoácido más grande)

io = PDBIO()
io.set_structure(estructura)
io.save("3BE3_mutado_TRP.pdb")