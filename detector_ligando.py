from Bio.PDB import PDBParser

def encontrar_ligandos_pdb(archivo_pdb):
    parser = PDBParser(QUIET=True)
    try:
        estructura = parser.get_structure('estructura_ligando', archivo_pdb)
        print(f"--- Buscando Ligandos en el Cristal: {archivo_pdb} ---")
        
        ligandos_encontrados = False
        
        # Recorremos la arquitectura del PDB
        for modelo in estructura:
            for cadena in modelo:
                for residuo in cadena:
                    # El id de un residuo en Biopython tiene 3 partes: (hetero_flag, res_number, insertion_code)
                    id_residuo = residuo.id
                    
                    # Si el hetero_flag empieza con 'H_', significa que es un LIGANDO (HETATM)
                    if id_residuo[0].startswith('H_'):
                        nombre_ligando = residuo.resname.strip()
                        numero_residuo = id_residuo[1]
                        
                        # Filtro: Ignoramos el agua convencional (HOH) para enfocarnos en ligandos reales
                        if nombre_ligando != "HOH":
                            ligandos_encontrados = True
                            print(f"\n¡LIGANDO DETECTADO!")
                            print(f"  Nombre Químico (Código): {nombre_ligando}")
                            print(f"  Ubicación en Cadena: [{cadena.id}]")
                            print(f"  Número de Residuo en PDB: {numero_residuo}")
                            print("-" * 50)
                            print("Copia y pega este comando en la consola de PyMOL para aislarlo:")
                            print(f"  select mi_ligando, resi {numero_residuo} and chain {cadena.id}")
                            print("  show sticks, mi_ligando")
                            print("  zoom mi_ligando, 8")
                            print("-" * 50)
                            
        if not ligandos_encontrados:
            print("No se detectaron ligandos (HETATM) de interés en este archivo PDB.")
            
    except FileNotFoundError:
        print(f"Error: No se encuentra el archivo '{archivo_pdb}'. Asegúrate de descargarlo.")
    except Exception as e:
        print(f"Ocurrió un error al analizar la estructura: {e}")

# --- CONFIGURACIÓN ---
# Descarga el archivo 1kya.pdb de rcsb.org y colócalo en la misma carpeta
ARCHIVO_PDB = "1KYA.pdb" 

if __name__ == "__main__":
    encontrar_ligandos_pdb(ARCHIVO_PDB)