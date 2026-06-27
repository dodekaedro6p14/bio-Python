from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

def localizar_bloque_perfecto(archivo_pdb, secuencia_del_html):
    # 1. Control de calidad de la entrada: Convertimos siempre a MAYÚSCULAS
    # y eliminamos espacios o guiones que se hayan colado al copiar del HTML
    secuencia_buscada = secuencia_del_html.upper().strip().replace("-", "")
    
    if not secuencia_buscada:
        print("Error: La secuencia de búsqueda está vacía.")
        return

    parser = PDBParser(QUIET=True)
    try:
        estructura = parser.get_structure('enzima', archivo_pdb)
        print(f"--- Escaneando Estructura PDB: {archivo_pdb} ---")
        print(f"Buscando el bloque patrón: '{secuencia_buscada}' (Procesado en Mayúsculas)\n")
        
        bloque_encontrado_global = False
        
        # 2. Recorremos los modelos y cadenas del archivo PDB
        for modelo in estructura:
            for cadena in modelo:
                secuencia_cadena = ""
                lista_residuos_num = []
                
                for residuo in cadena:
                    # Filtro: Solo tomamos aminoácidos estándar de la proteína
                    if residuo.id[0] == ' ':
                        # Convertimos a letra única (siempre devuelve mayúscula, ej: 'A')
                        letra = seq1(residuo.resname)
                        secuencia_cadena += letra
                        # Guardamos el número real identificador del residuo en el PDB
                        lista_residuos_num.append(residuo.id[1])
                
                # 3. Realizamos la búsqueda del bloque de texto en la cadena
                posicion_idx = secuencia_cadena.find(secuencia_buscada)
                
                if posicion_idx != -1:
                    bloque_encontrado_global = True
                    longitud = len(secuencia_buscada)
                    
                    # Extraemos los números de inicio y fin exactos del PDB
                    num_inicio = lista_residuos_num[posicion_idx]
                    num_fin = lista_residuos_num[posicion_idx + longitud - 1]
                    
                    print(f"¡ÉXITO! Bloque localizado en la CADENA: [{cadena.id}]")
                    print(f"Patrón detectado entre los residuos {num_inicio} y {num_fin}")
                    print("=" * 60)
                    print("Copia y pega esta línea en la consola de PyMOL:")
                    print(f"  select bloque_html, resi {num_inicio}-{num_fin}")
                    print("  color yellow, bloque_html")
                    print("  show sticks, bloque_html")
                    print("=" * 60 + "\n")
        
        if not bloque_encontrado_global:
            print("Resultado: No se encontró el bloque en ninguna cadena.")
            print("Verifica que el archivo PDB corresponda exactamente a la secuencia analizada.")
            
    except FileNotFoundError:
        print(f"Error: No se encontró el archivo '{archivo_pdb}' en esta carpeta.")
    except Exception as e:
        print(f"Ocurrió un error inesperado: {e}")

# --- CONFIGURACIÓN DE TU INVESTIGACIÓN ---
# Puedes usar letras minúsculas, mayúsculas o combinadas de tu HTML. El script lo resolverá solo.
ARCHIVO_PDB = "1Y2X.pdb"  
CADENA_HTML = "VteGdNlkaNLI"  # Ejemplo inicial real de la laccasa Trametes versicolor

if __name__ == "__main__":
    localizar_bloque_perfecto(ARCHIVO_PDB, CADENA_HTML)
