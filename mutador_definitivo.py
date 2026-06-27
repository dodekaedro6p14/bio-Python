from Bio.PDB import PDBParser, PDBIO

parser = PDBParser()
# Asegúrate de tener el archivo 3be3.pdb original en la misma carpeta
estructura = parser.get_structure('3BE3', '3be3.pdb')

# --- AQUÍ CAMBIAMOS LAS VARIABLES ---
POSICION = 30          # Cambiamos el 77 por el 30
NUEVO_AMINOACIDO = 'ASN' # Cambiamos 'TRP' por 'ASN' (Asparagina)
# ------------------------------------

# El script accede automáticamente a la cadena A y a la posición que le diste
residuo_a_mutar = estructura[0]['A'][POSICION]

# Hacemos una pequeña verificación de seguridad en la terminal
print(f"Residuo original detectado: {residuo_a_mutar.resname} en la posición {POSICION}")

# Aplicamos el cambio
residuo_a_mutar.resname = NUEVO_AMINOACIDO 
print(f"Mutación aplicada con éxito -> Nuevo nombre: {residuo_a_mutar.resname}")

# Guardamos con un nombre descriptivo para no confundirte en PyMOL
io = PDBIO()
io.set_structure(estructura)
io.save(f"3BE3_mutado_{NUEVO_AMINOACIDO}_{POSICION}.pdb")