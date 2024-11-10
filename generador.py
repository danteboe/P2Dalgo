import random
import string

def generar_peptidos(cantidad):
    """Genera una lista de péptidos aleatorios de 5 caracteres cada uno."""
    peptidos = []
    for _ in range(cantidad):
        peptido = ''.join(random.choices(string.ascii_uppercase, k=5))
        peptidos.append(peptido)
    return peptidos

def generar_celulas(n, max_coord, max_peptidos):
    """Genera una lista de células con atributos aleatorios."""
    celulas = []
    for i in range(1, n + 1):
        posX = random.randint(0, max_coord)
        posY = random.randint(0, max_coord)
        tipo = random.randint(1, 3)  # Tipo de célula: 1 (iniciadora), 2 (calculadora), 3 (ejecutora)
        peptidos = generar_peptidos(random.randint(1, max_peptidos))
        
        celula = f"{i} {posX} {posY} {tipo} " + " ".join(peptidos)
        celulas.append(celula)
    
    return celulas

def generar_caso_prueba(n, d, max_coord, max_peptidos):
    """Genera un caso de prueba completo con n células y distancia d."""
    caso = [f"{n} {d}"]
    celulas = generar_celulas(n, max_coord, max_peptidos)
    caso.extend(celulas)
    return "\n".join(caso)

def generar_pruebas(num_casos, n, d, max_coord, max_peptidos, output_file="pruebas_ultralisk.txt"):
    """Genera múltiples casos de prueba y los guarda en un archivo."""
    with open(output_file, "w") as f:
        f.write(f"{num_casos}\n")
        for _ in range(num_casos):
            caso_prueba = generar_caso_prueba(n, d, max_coord, max_peptidos)
            f.write(caso_prueba + "\n\n")  # Separador entre casos de prueba

# Parámetros de configuración
num_casos = 3          # Número de casos de prueba
n = 100000                 # Número de células por caso
d = 5                  # Distancia máxima para enviar mensajes
max_coord = 1000         # Máxima coordenada para las células
max_peptidos = 100000       # Máximo número de péptidos por célula

# Generar archivo de pruebas
generar_pruebas(num_casos, n, d, max_coord, max_peptidos)
print(f"Pruebas generadas en 'pruebas_ultralisk.txt'")
