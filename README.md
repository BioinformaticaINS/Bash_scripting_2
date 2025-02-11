# Bash scripting 2

## ¿Dónde se encuentra el archivo `.bashrc`?

El archivo `.bashrc` es un script de shell que se ejecuta cada vez que un usuario inicia una nueva sesión de terminal en sistemas operativos basados en Unix, como Linux y macOS. Este archivo se encuentra en el directorio **home** del usuario.

- **Ubicación**: `~/.bashrc`

### Explicación:
- `~` representa el directorio home del usuario actual. Por ejemplo, si el nombre de usuario es `usuario`, la ruta completa sería `/home/usuario/.bashrc` en sistemas Linux o `/Users/usuario/.bashrc` en macOS.

## ¿Cómo agregamos un directorio de forma PERMANENTE al PATH?

Para agregar un directorio de forma permanente al `PATH`, necesitas modificar el archivo `.bashrc` (o `.bash_profile` en algunos sistemas) y luego recargar la configuración.

## Pasos para agregar un directorio al PATH:

1. **Abrir el archivo `.bashrc` en un editor de texto**:
   ```bash
   nano ~/.bashrc
   ```
   (Puedes usar `nano`, `vim`, o cualquier otro editor de texto de tu preferencia).

2. **Agregar la siguiente línea al final del archivo**:
   ```bash
   export PATH="/ruta/al/directorio:$PATH"
   ```
   Reemplaza `"/ruta/al/directorio"` con la ruta del directorio que deseas agregar.

3. **Guardar y cerrar el archivo**:
   - En `nano`, presiona `CTRL + O` para guardar y `CTRL + X` para salir.

4. **Recargar el archivo `.bashrc` para aplicar los cambios**:
   ```bash
   source ~/.bashrc
   ```
   O simplemente cierra y vuelve a abrir la terminal.

### Ejemplo:

**Paso 1: Descarga del programa QUAST**

```bash
cd Bioprogrmas
wget https://downloads.sourceforge.net/project/quast/quast-5.3.0.tar.gz
tar xvfz quast-5.3.0.tar.gz
cd quast-5.3.0
ls
```


**Paso 2: Verificar la ubicación de QUAST**

Primero, asegúrate de que el ejecutable de QUAST esté en la carpeta correcta. Supongamos que has descargado y descomprimido QUAST en `~/Bioprograms/quast-5.3.0/`. El ejecutable de QUAST se encuentra en el directorio principal o en el subdirectorio `bin/`, dependiendo de la versión.

**Paso 2: Agregar QUAST al PATH** 

Para agregar QUAST al `PATH`, necesitas modificar el archivo `.bashrc` (o `.bash_profile` en algunos sistemas) y luego recargar la configuración.

**Pasos para agregar QUAST al PATH:**

1. **Abrir el archivo `.bashrc` en un editor de texto**:
   ```bash
   nano ~/.bashrc
   ```
   (Puedes usar `nano`, `vim`, o cualquier otro editor de texto de tu preferencia).

2. **Agregar la siguiente línea al final del archivo**:
   ```bash
   export PATH="$HOME/Bioprograms/quast-5.3.0:$PATH"
   ```
   Asegúrate de reemplazar `quast-5.3.0` con la versión exacta que tienes instalada si es diferente.

3. **Guardar y cerrar el archivo**:
   - En `nano`, presiona `CTRL + O` para guardar y `CTRL + X` para salir.

4. **Recargar el archivo `.bashrc` para aplicar los cambios**:
   ```bash
   source ~/.bashrc
   ```
   O simplemente cierra y vuelve a abrir la terminal.

### Paso 3: Verificar la instalación

Para asegurarte de que QUAST se ha agregado correctamente al `PATH`, puedes verificar la versión de QUAST ejecutando:

```bash
quast.py --version
```

¿Problemas? Vamos a corregirlos....

---

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
echo "Primer argumento (archivo FASTQ): $1" # (del 1 al 9) -> ${10} ${11}
echo "Número de argumentos: $#"
echo "PID del script: $$"
```

Este script podría usarse para procesar un archivo FASTQ, donde el primer argumento es el nombre del archivo.

---

## 2. **Comando `read`**

El comando `read` permite leer la entrada estándar del usuario y almacenarla en una variable.

Sintaxis: 

```
(base) ins_user@VirtualBox:~$ read var1 var 2 ...
```

**Ejemplo:**

```bash
(base) ins_user@VirtualBox:~$read var1
Hola # La palabra introducida se almacena en la
variable var1
(base) ins_user@VirtualBox:~$echo var1
Hola
(base) ins_user@VirtualBox:~$read var1 var2
analisis bioinformática
(base) ins_user@VirtualBox:~$echo var1
analisis
(base) ins_user@VirtualBox:~$echo var2
bioinformática
(base) ins_user@VirtualBox:~$read var1 var2
vamos incrementando nuestro conocimiento
(base) ins_user@VirtualBox:~$echo var1
vamos
(base) ins_user@VirtualBox:~$echo var2
incrementando nuestro conocimiento
```

**Ejemplo bioinformático:**
```bash
#!/bin/bash
echo "Introduce el nombre del archivo FASTA:"
read fasta_file
SRR123456
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
Pero ¿y si queremos buscar todos los números que aparecen en un texto? ¿O todas las líneas que empiezan por una letra mayúscula?


---

## 4. **Caracteres y Metacaracteres**

Los metacaracteres son caracteres especiales que tienen un significado específico en las expresiones regulares. Algunos de los más comunes son:

| **Metacarácter** | **Descripción**                                                                 | **Ejemplo**                                                                                     |
|-------------------|---------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------|
| `.`               | Coincide con cualquier carácter excepto el salto de línea (`\n`).               | `a.b` coincide con `"aab"`, `"acb"`, `"a1b"`, pero no con `"ab"` ni `"a\nb"`.                   |
| `*`               | Coincide con cero o más repeticiones del elemento anterior.                     | `ab*` coincide con `"a"`, `"ab"`, `"abb"`, `"abbb"`, etc.                                       |
| `+`               | Coincide con una o más repeticiones del elemento anterior.                      | `ab+` coincide con `"ab"`, `"abb"`, `"abbb"`, pero no con `"a"`.                                |
| `?`               | Coincide con cero o una repetición del elemento anterior.                       | `ab?` coincide con `"a"`, `"ab"`, pero no con `"abb"`.                                          |
| `[]`              | Define un conjunto de caracteres. Coincide con cualquiera de los caracteres dentro de los corchetes. | `[abc]` coincide con `"a"`, `"b"`, `"c"`, pero no con `"d"`.                                   |
| `[^]`             | Define un conjunto negado. Coincide con cualquier carácter que **no** esté en el conjunto. | `[^abc]` coincide con cualquier carácter excepto `"a"`, `"b"`, `"c"`.                          |
| `()`              | Agrupa patrones para aplicar operadores a un grupo completo.                    | `(ab)+` coincide con `"ab"`, `"abab"`, `"ababab"`, etc.                                         |
| `{}`              | Especifica un número exacto de repeticiones del elemento anterior.              | `a{3}` coincide con `"aaa"`. `a{2,4}` coincide con `"aa"`, `"aaa"`, `"aaaa"`.                   |
| `^`               | Coincide con el inicio de una línea o cadena.                                   | `^a` coincide con `"abc"`, pero no con `"bac"`.                                                 |
| `$`               | Coincide con el final de una línea o cadena.                                    | `a$` coincide con `"ba"`, pero no con `"abc"`.                                                  |
| `\`               | Escapa un metacarácter para tratarlo como un carácter literal.                  | `\.` coincide con un punto literal (`"."`), en lugar de interpretarlo como "cualquier carácter".|
| `\d`              | Coincide con cualquier dígito (equivalente a `[0-9]`).                          | `\d` coincide con `"1"`, `"2"`, etc., pero no con `"a"`.                                        |
| `\w`              | Coincide con cualquier carácter alfanumérico o guion bajo (equivalente a `[a-zA-Z0-9_]`). | `\w` coincide con `"a"`, `"1"`, `"_"`, pero no con `"@"`.                                      |
| `\s`              | Coincide con cualquier espacio en blanco (espacios, tabulaciones, saltos de línea). | `\s` coincide con `" "`, `"\t"`, `"\n"`.                                                       |
| `\D`, `\W`, `\S`  | Negación de `\d`, `\w`, `\s`. Coinciden con cualquier cosa que **no** sea un dígito, alfanumérico o espacio en blanco, respectivamente. | `\D` coincide con `"a"`, `"@"`, pero no con `"1"`.                                             |
| `|`               | Actúa como un operador OR lógico. Coincide con cualquiera de los patrones separados por `|`. | `a|b` coincide con `"a"` o `"b"`.                                                              |

### **Explicación Adicional**
1. **Agrupación con `()`**:
   - Los paréntesis permiten agrupar partes de una expresión regular para aplicar operadores como `*`, `+`, o `?` a todo el grupo.
   - Ejemplo: `(ab)+` significa que la secuencia `"ab"` debe repetirse una o más veces.

2. **Rangos en `[]`**:
   - Puedes usar rangos dentro de los corchetes para especificar conjuntos de caracteres.
   - Ejemplo: `[a-z]` coincide con cualquier letra minúscula, `[0-9]` coincide con cualquier dígito.

3. **Escapar metacaracteres con `\`**:
   - Si necesitas buscar un metacarácter como un carácter literal (por ejemplo, un punto `.`), debes escaparlo con una barra invertida (`\.`).

4. **Cuantificadores con `{}`**:
   - `{n}`: Exactamente `n` repeticiones.
   - `{n,}`: Al menos `n` repeticiones.
   - `{n,m}`: Entre `n` y `m` repeticiones.

**Ejemplo bioinformático:**

Creación de secuencias.fasta:


```
>Secuencia1
ATGCCGTAA
>Secuencia2
ATGGCTAGCTAA
>Secuencia3
ATGTTTAA
```

```bash
#!/bin/bash
# Buscar secuencias que comiencen con "ATG" y terminen con "TAA"
grep "^ATG.*TAA$" secuencias.fasta
```

---

## 5. **Rango de Caracteres**

Los rangos de caracteres permiten definir un conjunto de caracteres que pueden coincidir en una expresión regular.

### **Tabla de clases de caracteres**

| **Clase**       | **Caracteres**               | **Significado**                                                                 |
|------------------|------------------------------|---------------------------------------------------------------------------------|
| `[:alnum:]`      | `[A-Za-z0-9]`                | Caracteres alfanuméricos (letras mayúsculas, minúsculas y dígitos).             |
| `[:word:]`       | `[A-Za-z0-9_]`               | Caracteres alfanuméricos y el guion bajo (`_`).                                 |
| `[:alpha:]`      | `[A-Za-z]`                   | Caracteres alfabéticos (solo letras mayúsculas y minúsculas).                   |
| `[:blank:]`      | `[\t]`                       | Espacio y tabulador.                                                            |
| `[:space:]`      | `[\n\r\f\v]`                 | Espacios (incluye salto de línea, retorno de carro, etc.).                      |
| `[:digit:]`      | `[0-9]`                      | Dígitos (números del 0 al 9).                                                   |
| `[:lower:]`      | `[a-z]`                      | Letras minúsculas.                                                              |
| `[:upper:]`      | `[A-Z]`                      | Letras mayúsculas.                                                              |
| `[:punct:]`      | `[!@#$%^&*()_+=\-[\]{};':",./<>?~]` | Caracteres de puntuación.                                                       |

---

### **Tabla de rangos de caracteres**

| **Expresión** | **Coincidencia en el Texto**                                                                 |
|---------------|---------------------------------------------------------------------------------------------|
| `[abc]`       | El patrón coincide con la cadena si contiene una `a`, una `b` o una `c`.                    |
| `[a-c]`       | Equivalente a `[abc]`. Coincide con cualquier carácter entre `a` y `c`.                     |
| `c[aeo]sa`    | Coincide con las palabras `casa`, `cesa`, `cosa`.                                           |
| `[^abc]`      | El patrón coincide con la cadena si **no** contiene ninguna `a`, `b` o `c`.                 |
| `[0-9]`       | Coincide con una cadena que contenga cualquier dígito del `0` al `9`.                       |
| `[^0-9]`      | Coincide con una cadena que **no** contenga ningún dígito.                                  |

**Ejemplos bioinformáticos:**

### **Ejemplo 1: Uso de `[:alnum:]`**

**Objetivo**: Buscar líneas en un archivo FASTA que contengan solo caracteres alfanuméricos (letras y números).

**Archivo FASTA (`secuencias.fasta`)**:
```
>seq1
ATCGATCG
>seq2
12345678
>seq3
ATCG1234
>seq4
ATCG!@#$
```

**Comando**:
```bash
grep "^>[[:alnum:]]*$" secuencias.fasta
```

**Resultado**:
```
>seq1
>seq2
>seq3
```

**Explicación**: El comando busca líneas de encabezado que contienen solo letras y números. La línea `>seq4` no coincide porque contiene caracteres especiales (`!@#$`).

---

### **Ejemplo 2: Uso de `[:digit:]`**

**Objetivo**: Buscar secuencias que contengan números en su identificador.

**Archivo FASTA (`secuencias.fasta`)**:
```
>seq1
ATCGATCG
>seq2
12345678
>seq3
ATCG1234
>seq4
ATCG!@#$
```

**Comando**:
```bash
grep "[[:digit:]]" secuencias.fasta
```

**Resultado**:
```
>seq2
12345678
>seq3
ATCG1234
```

**Explicación**: El comando busca líneas que contienen dígitos. Las líneas `>seq2` y `>seq3` coinciden porque contienen números.

---

### **Ejemplo 3: Uso de `[A-Z]`**

**Objetivo**: Buscar secuencias que contengan solo letras mayúsculas (A, T, C, G).

**Archivo FASTA (`secuencias.fasta`)**:
```
>seq1
ATCGATCG
>seq2
12345678
>seq3
ATCG1234
>seq4
ATCG!@#$
```

**Comando**:
```bash
grep "^[A-Z]*$" secuencias.fasta
```

**Resultado**:
```
ATCGATCG
```

**Explicación**: El comando busca líneas que contienen solo letras mayúsculas. Solo la línea `ATCGATCG` coincide porque no contiene números ni caracteres especiales.

---

### **Ejemplo 4: Uso de `[^0-9]`**

**Objetivo**: Buscar secuencias que **no** contengan números.

**Archivo FASTA (`secuencias.fasta`)**:
```
>seq1
ATCGATCG
>seq2
12345678
>seq3
ATCG1234
>seq4
ATCG!@#$
```

**Comando**:
```bash
grep "^[^0-9]*$" secuencias.fasta
```

**Resultado**:
```
ATCGATCG
ATCG!@#$
```

**Explicación**: El comando busca líneas que no contienen números. Las líneas `ATCGATCG` y `ATCG!@#$` coinciden porque no tienen dígitos.

---

### **Ejemplo 5: Uso de `c[aeo]sa`**

**Objetivo**: Buscar nombres de genes que coincidan con el patrón `casa`, `cesa`, o `cosa`.

**Archivo de genes (`genes.txt`)**:
```
casa
cesa
cosa
cusa
cisa
```

**Comando**:
```bash
grep "c[aeo]sa" genes.txt
```

**Resultado**:
```
casa
cesa
cosa
```

**Explicación**: El comando busca nombres de genes que coinciden con el patrón `c[aeo]sa`. Los nombres `casa`, `cesa` y `cosa` coinciden, mientras que `cusa` y `cisa` no.

---

### **Ejemplo 6: Uso de `[:punct:]`**

**Objetivo**: Buscar líneas en un archivo de anotaciones que contengan caracteres de puntuación.

**Archivo de anotaciones (`anotaciones.txt`)**:
```
Gene1: Expressed in liver.
Gene2: Expressed in brain, heart.
Gene3: Expressed in kidney.
Gene4: Expressed in lung!
```

**Comando**:
```bash
grep "[[:punct:]]" anotaciones.txt
```

**Resultado**:
```
Gene1: Expressed in liver.
Gene2: Expressed in brain, heart.
Gene4: Expressed in lung!
```

**Explicación**: El comando busca líneas que contienen caracteres de puntuación. Las líneas `Gene1`, `Gene2` y `Gene4` coinciden porque contienen `:`, `,` y `!`.

---

### **Ejemplo 7: Uso de `[0-9]`**

**Objetivo**: Buscar secuencias que contengan dígitos.

**Archivo FASTA (`secuencias.fasta`)**:
```
>seq1
ATCGATCG
>seq2
12345678
>seq3
ATCG1234
>seq4
ATCG!@#$
```

**Comando**:
```bash
grep "[0-9]" secuencias.fasta
```

**Resultado**:
```
12345678
ATCG1234
```

**Explicación**: El comando busca líneas que contienen dígitos. Las líneas `12345678` y `ATCG1234` coinciden porque contienen números.

---

### **Ejemplo 8: Uso de `[:upper:]`**

**Objetivo**: Buscar secuencias que contengan solo letras mayúsculas (A, T, C, G).

**Archivo FASTA (`secuencias.fasta`)**:
```
>seq1
ATCGATCG
>seq2
12345678
>seq3
ATCG1234
>seq4
ATCG!@#$
```

**Comando**:
```bash
grep "^[[:upper:]]*$" secuencias.fasta
```

**Resultado**:
```
ATCGATCG
```

**Explicación**: El comando busca líneas que contienen solo letras mayúsculas. Solo la línea `ATCGATCG` coincide porque no contiene números ni caracteres especiales.

Este comando busca secuencias que contienen solo los nucleótidos A, T, C y G.

---

## 6. **Cuantificadores**

Los cuantificadores indican cuántas veces debe aparecer un carácter o grupo de caracteres.

- `?`: Cero o una vez.
- `*`: Cero o más veces.
- `+`: Una o más veces.
- `{n}`: Exactamente `n` veces.
- `{n,m}`: Entre `n` y `m` veces.

¡Claro! Vamos a desglosar el uso del comando `grep "ATG\{3,5\}" secuencias.fasta` con un ejemplo concreto y resultados esperados.

---

### **Explicación del Comando**

El comando `grep "ATG\{3,5\}" secuencias.fasta` busca en el archivo `secuencias.fasta` patrones que contengan la secuencia "ATG" repetida entre **3 y 5 veces**. Aquí está el desglose:

- `ATG`: Es la secuencia que estamos buscando.
- `\{3,5\}`: Es un cuantificador que indica que "ATG" debe aparecer entre 3 y 5 veces consecutivas.

---

### **Ejemplo con Datos Biológicos**
### **Aplicación Bioinformática**

Este tipo de búsqueda es útil en bioinformática para identificar regiones específicas en secuencias de ADN. Por ejemplo:
- **Codones de inicio**: "ATG" es el codón de inicio en el código genético. Buscar repeticiones de "ATG" podría ayudar a identificar regiones de interés en secuencias genómicas.
- **Motivos repetidos**: Algunas secuencias repetidas pueden tener funciones biológicas específicas, como la regulación de la expresión génica.

---

### **Ejemplo Ampliado**

Supongamos que tienes un archivo FASTA más grande y quieres buscar repeticiones de "ATG" en las secuencias:

#### **Archivo FASTA (`secuencias.fasta`)**:
```
>seq1
ATGATGATGCGTAGCTAG
>seq2
ATGATGATGATGATGCGTAG
>seq3
ATGATGATGATGATGATG
>seq4
ATGCGTAGCTAG
>seq5
ATGATG
>seq6
ATGATGATG
```

#### **Comando**:
```bash
grep "ATG\{3,5\}" secuencias.fasta
```

#### **Resultado Esperado**:
```
ATGATGATGCGTAGCTAG
ATGATGATGATGATGCGTAG
ATGATGATGATGATGATG
ATGATGATG
```
---

## 7. **Listas con BASH**

Las listas (arrays) en Bash permiten almacenar múltiples valores en una sola variable.

• Variable con múltiples valores (distinta naturaleza).
• Valores secuenciales.
• NO hay límite máximo.
• El índice del campo comienza en cero.

**Ejemplo:**

```bash
(base) ins_user@VirtualBox:~$numeros=(1 2 3 4 5)
(base) ins_user@VirtualBox:~$num[0]=1
(base) ins_user@VirtualBox:~$num[1]=2
(base) ins_user@VirtualBox:~$num[2]=3
```

Mostrar el contenido completo de una lista

Sintaxis: `echo ${lista[@]}`
        : `echo ${lista[*]}`


```bash
(base) ins_user@VirtualBox:~$numeros=(1 2 3 4 5)
(base) ins_user@VirtualBox:~$echo ${numeros[*]}
1 2 3 4 5
(base) ins_user@VirtualBox:~$echo ${numeros[@]}
```

Acceder a un elemento de una lista

sintaxis: `echo ${numeros[indice o posicion]}`

```bash
(base) ins_user@VirtualBox:~$numeros=(1 2 3 4 5)
(base) ins_user@VirtualBox:~${numeros[0]}
1
(base) ins_user@VirtualBox:~${numeros[1]}
2
....
```

Añadir o modificar un elemento a una lista

```bash
(base) ins_user@VirtualBox:~$numeros=(1 2 3 4 5)
(base) ins_user@VirtualBox:~${numeros[*]}
1 2 3 4 5
(base) ins_user@VirtualBox:~$numeros[5]=6
(base) ins_user@VirtualBox:~${numeros[*]}
1 2 3 4 5 6
```

Añadir al final, sin necesidad de especificar el índice

Sintaxis: `numeros+=(nuevo_valor(es))`

```bash
(base) ins_user@VirtualBox:~$numeros=(1 2 3 4 5)
(base) ins_user@VirtualBox:~${numeros[*]}
1 2 3 4 5
(base) ins_user@VirtualBox:~$numeros+=(8 9 10)
(base) ins_user@VirtualBox:~${numeros[*]}
1 2 3 4 5 8 9 10
```

Longitud de un lista 

sintaxis: `echo ${#numeros[*]}`

```bash
(base) ins_user@VirtualBox:~$numeros=(1 2 3 4 5)
(base) ins_user@VirtualBox:~${numeros[*]}
1 2 3 4 5
(base) ins_user@VirtualBox:~$echo ${#numeros[*]}
5
```

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

![](https://procomsys.wordpress.com/wp-content/uploads/2018/05/escon.png)

![](https://bookdown.org/hneth/ds4psy/images/bart/bart_board_for.png)

**Sintaxis de bucle FOR**

Se usa para ejecutar un conjunto de comandos dado un número conocido de veces.

```bash
for i in 1 2 3 4 5 .. N # i iterador, N: Lista de elementos 
do
<-->COMANDO 1 # <--> sangrado
<-->COMANDO 2
    ...
<-->COMANDO N
done
```

Ejemplos:

```bash
#!/bin/bash

for num in 1 2 3 4 5
do
   echo "Hello $num"
done
```

Por línea de comandos

```bash
(base) ins_user@VirtualBox:~$ for num in 1 2 3 4 5; do echo "Hello $num"; done
```

Aplicado a lista de números:

```bash
#!/bin/bash
for num in {1..5}
do
   echo "Hello $num"
done
```

Por línea de comandos

```bash
(base) ins_user@VirtualBox:~$ for num in {1..10..2}; do echo "$num"; done
```

```bash
#!/bin/bash
num=(1 2 3 4 5)
for num in ${num[@]}
do
   echo "Hello $num"
done
```

```bash
#!/bin/bash

genes=(SOX13 PAX5 TC1 ADF)
for gene in ${genes[@]}
do
   echo "La longituf de $gene es ${#gene}"
done
```

## Ejemplo de Aplicación Bioinformática

```bash
#!/bin/bash
# Script para automatizar el procesamiento de archivos FASTQ con datos de Ebola

# Definir variables para las carpetas
REF_DIR="Proyecto_NGS/refs"          # Carpeta para el genoma de referencia
RAW_DATA_DIR="Proyecto_NGS/raw_data" # Carpeta para los datos crudos (FASTQ)
RESULTS_DIR="Proyecto_NGS/results"   # Carpeta para los resultados (BAM, SAM, etc.)

# Crear directorios si no existen
mkdir -p $REF_DIR
mkdir -p $RAW_DATA_DIR
mkdir -p $RESULTS_DIR

# Instalar herramientas necesarias
echo "Instalando herramientas..."
sudo apt update && sudo apt upgrade
sudo apt install -y bwa bowtie2 samtools fastp
sudo pip install bio --upgrade

# Descargar el genoma de referencia (Ebola, cepa de 1976)
echo "Descargando genoma de referencia..."
bio fetch AF086833 --format fasta > $REF_DIR/ebola_ref.fa

# Descargar datos de secuenciación (SRR1972739, 10,000 lecturas)
echo "Descargando datos de secuenciación..."
fastq-dump -X 10000 --split-files SRR1972739 -O $RAW_DATA_DIR

# Crear índice con BWA
echo "Creando índice con BWA..."
bwa index $REF_DIR/ebola_ref.fa

# Crear índice con Bowtie2
echo "Creando índice con Bowtie2..."
bowtie2-build $REF_DIR/ebola_ref.fa $REF_DIR/ebola_ref

# Verificar los archivos generados
echo "Verificando archivos generados..."
ls $REF_DIR/

# Lista de archivos FASTQ (lecturas forward y reverse)
fastq_files=("$RAW_DATA_DIR/SRR1972739_1.fastq" "$RAW_DATA_DIR/SRR1972739_2.fastq")

# Procesar cada archivo FASTQ con fastp
echo "Procesando archivos FASTQ con fastp..."
for file in "${fastq_files[@]}"; do
    echo "Procesando $file..."
    
    # Paso 1: Control de calidad con FastQC
    echo "Realizando control de calidad con FastQC..."
    fastqc $file -o $RESULTS_DIR/fastqc_output
    
    # Paso 2: Filtrar secuencias de baja calidad con fastp
    echo "Filtrando secuencias de baja calidad con fastp..."
    fastp -i $file \
        -o "${file%.fastq}_filtered.fastq" \
        --html "${file%.fastq}_fastp_report.html" \
        --json "${file%.fastq}_fastp_report.json" \
        --qualified_quality_phred 20 \
        --length_required 50
    
    echo "Filtrado de $file completado."
done

# Paso 3: Alineamiento con BWA (modo paired-end)
echo "Alineando lecturas con BWA..."
bwa mem $REF_DIR/ebola_ref.fa \
    $RAW_DATA_DIR/SRR1972739_1_filtered.fastq \
    $RAW_DATA_DIR/SRR1972739_2_filtered.fastq > $RESULTS_DIR/bwa_output.sam

# Verificar resultados iniciales
echo "Verificando resultados de alineación..."
head -n 20 $RESULTS_DIR/bwa_output.sam

# Convertir SAM a BAM y ordenar
echo "Convirtiendo SAM a BAM y ordenando..."
samtools view -S -b $RESULTS_DIR/bwa_output.sam > $RESULTS_DIR/bwa_output.bam
samtools sort $RESULTS_DIR/bwa_output.bam -o $RESULTS_DIR/bwa_output_sorted.bam

# Indexar el archivo BAM para visualización en IGV
echo "Indexando archivo BAM..."
samtools index $RESULTS_DIR/bwa_output_sorted.bam

echo "Procesamiento completado. Archivos BAM listos para visualización en IGV."
echo "Resultados guardados en la carpeta: $RESULTS_DIR"
```

---

### **Explicación de las Variables**

1. **`REF_DIR="refs"`**:
   - Define la carpeta donde se almacenará el genoma de referencia y sus índices.
   - Ejemplo: `refs/ebola_ref.fa`, `refs/ebola_ref.1.bt2`, etc.

2. **`RAW_DATA_DIR="raw_data"`**:
   - Define la carpeta donde se descargarán y almacenarán los archivos FASTQ crudos.
   - Ejemplo: `raw_data/SRR1972739_1.fastq`, `raw_data/SRR1972739_2.fastq`.

3. **`RESULTS_DIR="results"`**:
   - Define la carpeta donde se guardarán los resultados del procesamiento, como archivos BAM, SAM, y salidas de FastQC.
   - Ejemplo: `results/bwa_output.sam`, `results/bwa_output_sorted.bam`, `results/fastqc_output/`.

---

### **Estructura de Carpetas Resultante**

Después de ejecutar el script, la estructura de carpetas será la siguiente:

```
.
├── refs/
│   ├── ebola_ref.fa
│   ├── ebola_ref.fa.amb
│   ├── ebola_ref.fa.ann
│   ├── ebola_ref.fa.bwt
│   ├── ebola_ref.fa.pac
│   ├── ebola_ref.fa.sa
│   ├── ebola_ref.1.bt2
│   ├── ebola_ref.2.bt2
│   ├── ebola_ref.3.bt2
│   ├── ebola_ref.4.bt2
│   ├── ebola_ref.rev.1.bt2
│   └── ebola_ref.rev.2.bt2
├── raw_data/
│   ├── SRR1972739_1.fastq
│   └── SRR1972739_2.fastq
└── results/
    ├── fastqc_output/
    │   ├── SRR1972739_1_fastqc.html
    │   └── SRR1972739_2_fastqc.html
    ├── SRR1972739_1_filtered.fastq
    ├── SRR1972739_2_filtered.fastq
    ├── SRR1972739_1_fastp_report.html
    ├── SRR1972739_1_fastp_report.json
    ├── SRR1972739_2_fastp_report.html
    ├── SRR1972739_2_fastp_report.json
    ├── bwa_output.sam
    ├── bwa_output.bam
    ├── bwa_output_sorted.bam
    └── bwa_output_sorted.bam.bai
```

![MUCHAS GRACIAS](https://pbs.twimg.com/media/BdGmyPqCQAEBIeX.jpg)
