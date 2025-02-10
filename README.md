# Bash scripting 2

## 1. **Variables Especiales**

Las variables especiales en Bash son predefinidas y tienen significados específicos. Algunas de las más comunes son:

- `$0`: Nombre del script.
- `$1`, `$2`, ..., `$9`: Argumentos pasados al script.
- `$#`: Número de argumentos pasados al script.
- `$$`: PID (Process ID) del script actual.
- `$?`: Estado de salida del último comando ejecutado.

**Ejemplo bioinformático:**
```bash
#!/bin/bash
echo "Nombre del script: $0"
echo "Primer argumento (archivo FASTQ): $1"
echo "Número de argumentos: $#"
echo "PID del script: $$"
```

Este script podría usarse para procesar un archivo FASTQ, donde el primer argumento es el nombre del archivo.

---

## 2. **Comando `read`**

El comando `read` permite leer la entrada del usuario y almacenarla en una variable.

**Ejemplo bioinformático:**
```bash
#!/bin/bash
echo "Introduce el nombre del archivo FASTA:"
read fasta_file
echo "Procesando el archivo $fasta_file..."
```

Este script solicita al usuario el nombre de un archivo FASTA y luego lo procesa.

---

## 3. **Expresiones Regulares**

Las expresiones regulares (regex) son patrones que permiten buscar y manipular texto. Son útiles para filtrar y procesar datos en archivos de texto.

**Ejemplo bioinformático:**
```bash
#!/bin/bash
# Buscar todas las secuencias que contengan "ATG" en un archivo FASTA
grep "ATG" secuencias.fasta
```

Este comando busca todas las líneas en un archivo FASTA que contengan el codón de inicio "ATG".

---

## 4. **Caracteres y Metacaracteres**

Los metacaracteres son caracteres especiales que tienen un significado específico en las expresiones regulares. Algunos de los más comunes son:

- `.`: Cualquier carácter.
- `*`: Cero o más repeticiones del carácter anterior.
- `+`: Una o más repeticiones del carácter anterior.
- `[]`: Define un rango de caracteres.

**Ejemplo bioinformático:**
```bash
#!/bin/bash
# Buscar secuencias que comiencen con "ATG" y terminen con "TAA"
grep "^ATG.*TAA$" secuencias.fasta
```

Este comando busca secuencias que comienzan con "ATG" y terminan con "TAA".

---

## 5. **Rango de Caracteres**

Los rangos de caracteres permiten definir un conjunto de caracteres que pueden coincidir en una expresión regular.

**Ejemplo bioinformático:**
```bash
#!/bin/bash
# Buscar secuencias que contengan solo nucleótidos (A, T, C, G)
grep "^[ATCG]*$" secuencias.fasta
```

Este comando busca secuencias que contienen solo los nucleótidos A, T, C y G.

---

## 6. **Cuantificadores**

Los cuantificadores indican cuántas veces debe aparecer un carácter o grupo de caracteres.

- `?`: Cero o una vez.
- `*`: Cero o más veces.
- `+`: Una o más veces.
- `{n}`: Exactamente `n` veces.
- `{n,m}`: Entre `n` y `m` veces.

**Ejemplo bioinformático:**
```bash
#!/bin/bash
# Buscar secuencias que contengan entre 3 y 5 repeticiones de "ATG"
grep "ATG\{3,5\}" secuencias.fasta
```

Este comando busca secuencias que contienen entre 3 y 5 repeticiones del codón "ATG".

---

## 7. **Partes Básicas de un Shell Script**

Un script de Shell tiene una estructura básica que incluye:

- **Shebang (`#!/bin/bash`)**: Indica que el script debe ser ejecutado con Bash.
- **Comentarios**: Líneas que comienzan con `#` y son ignoradas por el intérprete.
- **Código**: Comandos y lógica del script.

**Ejemplo bioinformático:**
```bash
#!/bin/bash
# Script para contar el número de secuencias en un archivo FASTA
echo "Número de secuencias en el archivo:"
grep -c "^>" secuencias.fasta
```

Este script cuenta el número de secuencias en un archivo FASTA.

---

## 8. **Listas con BASH**

Las listas (arrays) en Bash permiten almacenar múltiples valores en una sola variable.

**Ejemplo bioinformático:**
```bash
#!/bin/bash
# Crear una lista de nombres de archivos FASTQ
fastq_files=("sample1.fastq" "sample2.fastq" "sample3.fastq")

# Procesar cada archivo FASTQ
for file in "${fastq_files[@]}"; do
    echo "Procesando $file..."
    # Aquí iría el comando para procesar el archivo FASTQ
done
```

Este script procesa una lista de archivos FASTQ.

---

## 9. **Control de Flujo**

El control de flujo permite ejecutar comandos condicionalmente o en bucles.

**Ejemplo bioinformático:**
```bash
#!/bin/bash
# Bucle para procesar múltiples archivos FASTQ
for file in *.fastq; do
    if [[ -f "$file" ]]; then
        echo "Procesando $file..."
        # Aquí iría el comando para procesar el archivo FASTQ
    else
        echo "$file no es un archivo válido."
    fi
done
```

Este script procesa todos los archivos FASTQ en el directorio actual.

---

## Ejemplo de Aplicación Bioinformática

**Procesamiento Automatizado de Archivos FASTQ:**

```bash
#!/bin/bash
# Script para automatizar el procesamiento de archivos FASTQ

# Lista de archivos FASTQ
fastq_files=("sample1.fastq" "sample2.fastq" "sample3.fastq")

# Procesar cada archivo FASTQ
for file in "${fastq_files[@]}"; do
    echo "Procesando $file..."
    
    # Paso 1: Control de calidad con FastQC
    fastqc $file
    
    # Paso 2: Filtrar secuencias de baja calidad
    trimmomatic SE -phred33 $file "${file%.fastq}_filtered.fastq" LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50
    
    # Paso 3: Alineamiento con BWA
    bwa mem referencia.fasta "${file%.fastq}_filtered.fastq" > "${file%.fastq}.sam"
    
    echo "Procesamiento de $file completado."
done
```

Este script automatiza el procesamiento de archivos FASTQ, incluyendo control de calidad, filtrado y alineamiento.

---

### Conclusión

El dominio de estos conceptos de Bash scripting es esencial para la automatización de tareas bioinformáticas, como el procesamiento de archivos FASTQ y FASTA. La combinación de expresiones regulares, listas y control de flujo permite crear scripts potentes y eficientes para el análisis de datos biológicos.
