
# GUÍA DE USO Y FLUJO DE TRABAJO EN EL CLÚSTER (VirResp)

---

## 🎯 EJEMPLO DE FLUJO DE TRABAJO (LO QUE VAMOS A HACER)

Para aprender a usar el sistema rápidamente, realizaremos un análisis filogenético completo a partir de un archivo de secuencias llamado `dummy.fasta`.

El flujo de trabajo estándar en el clúster consta de 3 fases principales:


```text
[Secuencias raw] ──> 1. CD-HIT ──> 2. MAFFT ──> 3. IQ-TREE ──> [Árbol Filogenético]
                  (Limpieza)     (Alineamiento)   (Reconstrucción)
```


1. **CD-HIT**: Elimina secuencias duplicadas o altamente idénticas para reducir redundancia.
2. **MAFFT**: Alinea las secuencias de nucleótidos o aminoácidos.
3. **IQ-TREE**: Construye el árbol filogenético por Máxima Verosimilitud (Maximum Likelihood).

*(Nota: En la carpeta de plantillas del repositorio en el clúster `/data/cnm/vrg/VirResp-cluster` ya están creados y optimizados todos los scripts. **No tienes que modificar ningún archivo ni script interno**, únicamente copiarlos a tu carpeta de trabajo y ejecutarlos).*

---

## ⚡ GUÍA DE USO RÁPIDO (PASO A PASO)

Sigue esta secuencia de comandos de inicio a fin. Copia y pega directamente ajustando solo tu usuario.

---

### 🔷 PASO 1: LOCALIZAR TUS ARCHIVOS EN LOCAL

> 🟦 **ENTORNO: TU ORDENADOR (LOCAL)**

Abre la terminal de tu ordenador y navega hasta la carpeta donde tienes tu archivo de trabajo (por ejemplo `dummy.fasta`).


```bash
# Navegar a la carpeta donde tienes tu archivo
cd /ruta/a/tu/carpeta/local

# Comprobar que el archivo dummy.fasta está ahí
ls
```


---

### 🔷 PASO 2: SUBIR ARCHIVOS AL CLÚSTER CON RSYNC

> 🟦 **ENTORNO: TU ORDENADOR (LOCAL)**

Envía tu archivo `dummy.fasta` (o `.xml`) desde tu ordenador al clúster. Al ejecutar el comando te pedirá tu contraseña de usuario del clúster.


```bash
# Sintaxis general:
# rsync -avz -e "ssh -p 32122" ./archivo.extension usuario@portutatis.isciii.es:/data/cnm/vrg/usuario/

# Ejemplo práctico enviando el archivo dummy.fasta:
rsync -avz -e "ssh -p 32122" ./dummy.fasta usuario@portutatis.isciii.es:/data/cnm/vrg/usuario/
```


---

### 🟩 PASO 3: ABRIR TERMINAL Y ACCEDER AL CLÚSTER

> 🟩 **ENTORNO: CLÚSTER (PORTUTATIS)**

Abre una nueva ventana o pestaña de la terminal para conectarte al clúster por SSH e introduce tu contraseña.


```bash
ssh -p 32122 usuario@portutatis.isciii.es
```


---

### 🟩 PASO 4: NAVEGAR A TU CARPETA Y VERIFICAR EL ARCHIVO

> 🟩 **ENTORNO: CLÚSTER (PORTUTATIS)**

Una vez dentro del clúster, muévete a tu directorio personal de trabajo y comprueba que `dummy.fasta` ha llegado correctamente.


```bash
# Ir a tu carpeta de trabajo en el clúster
cd /data/cnm/vrg/usuario/

# Listar los archivos para confirmar que dummy.fasta está dentro
ls -l
```


---

### 🟩 PASO 5: EJECUTAR EL FLUJO DE TRABAJO (CD-HIT ➔ MAFFT ➔ IQ-TREE)

> 🟩 **ENTORNO: CLÚSTER (PORTUTATIS)**

Los scripts maestro y plantillas están preconfigurados en la ruta `/data/cnm/vrg/VirResp-cluster`. No cambies nada de esos archivos. Para lanzar cada tarea, simplemente ejecuta el comando máster correspondiente pasando como parámetro tu archivo de entrada.

#### 5.1. CD-HIT (Filtrado de redundancia)
Lanza el filtrado de secuencias sobre `dummy.fasta`. Generará el archivo de salida `dummy_cdhit.fasta`.


```bash
/data/cnm/vrg/VirResp-cluster/scripts/master_cdhit.sh dummy.fasta
```

> 🔗 [*Ver explicación detallada y parámetros de la plantilla CD-HIT*](#-plantilla-cd-hit)

#### 5.2. MAFFT (Alineamiento de secuencias)
Una vez terminado CD-HIT, lanza el alineamiento sobre el archivo filtrado `dummy_cdhit.fasta`. Generará `dummy_cdhit_mafft.fasta`.


```bash
/data/cnm/vrg/VirResp-cluster/scripts/master_mafft.sh dummy_cdhit.fasta
```

> 🔗 [*Ver explicación detallada y parámetros de la plantilla MAFFT*](#-plantilla-mafft)

#### 5.3. IQ-TREE (Reconstrucción Filogenética)
Lanza IQ-TREE sobre el alineamiento `dummy_cdhit_mafft.fasta`. Generará el árbol filogenético y archivos asociados.


```bash
/data/cnm/vrg/VirResp-cluster/scripts/master_iqtree.sh dummy_cdhit_mafft.fasta
```

> 🔗 [*Ver explicación detallada y parámetros de la plantilla IQ-TREE*](#-plantilla-iq-tree)

#### 5.4. BEAST (Análisis Bayesiano - Opcional)
Si lo que deseas ejecutar es un análisis Bayesiano con BEAST a partir de un archivo `.xml` preparado previamente con BEAUti:


```bash
/data/cnm/vrg/VirResp-cluster/scripts/master_beast.sh mi_analisis.xml
```

> 🔗 [*Ver explicación detallada y uso de BEAST / BEAUti*](#-plantilla-beast-y-uso-de-beauti)

---

### 🟩 PASO 6: CONSULTAR EL ESTADO DE TUS TRABAJOS EN SLURM

> 🟩 **ENTORNO: CLÚSTER (PORTUTATIS)**

Mientras las tareas se ejecutan en segundo plano, puedes comprobar su progreso con estos sencillos comandos:


```bash
# Ver el estado de tus tareas en la cola de procesamiento
squeue -u tu_usuario

# Ver detalles o consumo de recursos de un trabajo específico (reemplaza ID_TRABAJO)
scontrol show job ID_TRABAJO

# Cancelar o detener una tarea si te has equivocado
scancel ID_TRABAJO
```


#### ¿Cómo interpretar los estados en `squeue`?
* **`R` (Running)**: La tarea se está ejecutando correctamente.
* **`PD` (Pending)**: La tarea está en espera en la cola hasta que haya nodos/recursos libres.
* **`CG` (Completing)**: La tarea está finalizando y guardando los resultados.

---

### 🔷 PASO 7: DESCARGAR LOS RESULTADOS A TU ORDENADOR

> 🟦 **ENTORNO: TU ORDENADOR (LOCAL)**

Una vez finalizados los análisis en el clúster, vuelve a la terminal de tu ordenador (LOCAL) para traer todos los resultados generados de vuelta a tu equipo.


```bash
# Sintaxis general:
# rsync -avz -e "ssh -p 32122" usuario@portutatis.isciii.es:/data/cnm/vrg/usuario/ ./carpeta_destino_local/

# Ejemplo para traer toda tu carpeta de trabajo al directorio actual en local:
rsync -avz -e "ssh -p 32122" usuario@portutatis.isciii.es:/data/cnm/vrg/usuario/* ./resultados_cluster/
```


---

## 🛠️ COMANDOS BÁSICOS DE LINUX (IMPRESCINDIBLES)

| Comando | Función / Uso |
| :--- | :--- |
| `cd /ruta/` | Cambiar de directorio / carpeta. |
| `cd ..` | Subir un nivel de carpeta hacia atrás. |
| `ls` | Listar los archivos y carpetas de la ubicación actual. |
| `ls -lh` | Listar archivos mostrando su tamaño en formato legible (MB, GB). |
| `pwd` | Muestra la ruta completa del directorio donde estás situado. |
| `mkdir nombre_carpeta` | Crear una carpeta nueva. |
| `cp origen destino` | Copiar un archivo a otra ubicación. |
| `mv origen destino` | Mover o cambiar de nombre un archivo o carpeta. |
| `cat archivo.txt` | Muestra todo el contenido del archivo en la pantalla. |
| `head -n 20 archivo.fasta` | Muestra solo las primeras 20 líneas del archivo. |

---

### ⚠️ AVISO IMPORTANTE SOBRE ELIMINACIÓN DE ARCHIVOS

> 🛑 **¡ATENCIÓN! EN EL CLÚSTER NO EXISTE PAPELERA DE RECICLAJE**
> 
> * **`rm archivo`**: Elimina un archivo permanentemente.
> * **`rm -r carpeta`**: Elimina una carpeta y todo su contenido de forma recursiva.
> 
> **NO ES POSIBLE RECUPERAR NINGÚN ARCHIVO O CARPETA UNA VEZ BORRADO.**
> Revisa minuciosamente la ruta y el nombre del archivo antes de presionar `Enter` al usar `rm` o `rm -r`.

---

## 📄 EXPLICACIÓN DE PLANTILLAS Y HERRAMIENTAS

A continuación se detalla el funcionamiento interno de cada comando para virólogos que requieran entender qué ejecuta la plantilla por debajo.

---

### 🧬 PLANTILLA CD-HIT

Reduce la redundancia de secuencias mediante agrupamiento (clustering) por similitud de identidad.


```bash
# Comando interno ejecutado en el script lanzar_cdhit.sh:
cd-hit -i entrada.fasta -o salida_cdhit.fasta -c 0.99 -n 5 -M 16000 -T 4
```


* **`-i`**: Archivo FASTA de entrada.
* **`-o`**: Nombre del archivo FASTA de salida filtrado.
* **`-c 0.99`**: Umbral de identidad de secuencia (99%). Secuencias con una similitud de nucleótidos $\ge 99\%$ se agrupan en una sola representación.
* **`-n 5`**: Tamaño de palabra (word size). Recomendado 5 para el rango de identidad $0.90 - 1.00$.
* **`-M 16000`**: Límite de memoria asignada en MB (16 GB).
* **`-T 4`**: Número de hilos de procesamiento (threads) asignados.

---

### 🧬 PLANTILLA MAFFT

Alineamiento múltiple de secuencias de nucleótidos o proteínas de gran precisión.


```bash
# Comando interno ejecutado en el script lanzar_mafft.sh:
mafft --auto --thread 8 entrada_cdhit.fasta > salida_aligned.fasta
```


* **`--auto`**: Selecciona automáticamente la estrategia de alineamiento óptima (L-INS-i, FFT-NS-i, etc.) en función del número y la longitud de las secuencias.
* **`--thread 8`**: Utiliza 8 hilos de procesamiento en el nodo asignado para acelerar el cálculo.
* **`> salida_aligned.fasta`**: Redirige la salida del alineamiento hacia el archivo indicado.

---

### 🧬 PLANTILLA IQ-TREE

Construcción de árboles filogenéticos por Máxima Verosimilitud (Maximum Likelihood) con selección automática del modelo de sustitución.


```bash
# Comando interno ejecutado en el script lanzar_iqtree.sh:
iqtree -s salida_aligned.fasta -m MFP -bb 1000 -nt AUTO
```


* **`-s`**: Archivo de alineamiento de entrada (formato FASTA, PHYLIP o NEXUS).
* **`-m MFP`**: *ModelFinder Plus*. Evalúa automáticamente todos los modelos de sustitución nucleotídica y elige el de mejor ajuste estadístico.
* **`-bb 1000`**: Realiza 1000 réplicas de *Ultrafast Bootstrap* (UFBoot) para evaluar el soporte de ramas y nodos del árbol.
* **`-nt AUTO`**: Determina automáticamente el número óptimo de hilos de CPU en el nodo.

---

### 🧬 PLANTILLA BEAST Y USO DE BEAUTI

BEAST es una herramienta de inferencia filogenética Bayesiana orientada a estimar fechas de divergencia, tasas evolutivas y dinámicas poblacionales a lo largo del tiempo.


```bash
# Comando interno ejecutado en el script lanzar_beast.sh:
beast -threads 8 mi_analisis.xml
```


#### Flujo de trabajo completo para análisis con BEAST


```text
[BEAUti] ──> Genera .xml ──> [BEAST en Clúster] ──> Genera .log y .trees ──> [Tracer / TreeAnnotator / FigTree]
```


1. **Creación del XML (BEAUti)**: Abre BEAUti en tu ordenador local (o en el entorno configurado), carga tu alineamiento, define las fechas de muestreo, el modelo de sustitución, el reloj molecular y el modelo poblacional, y guarda el archivo resultante como `analisis.xml`.
2. **Ejecución de la cadena (BEAST)**: Sube `analisis.xml` al clúster y lanza la plantilla de BEAST. El clúster procesará las cadenas MCMC y generará dos archivos clave de salida:
   * `analisis.log`: Contiene la traza numérica de los parámetros muestreados.
   * `analisis.trees`: Contiene la colección de árboles generados durante la simulación.
3. **Evaluación de parámetros (Tracer)**: Una vez descargado el archivo `analisis.log` a tu ordenador, abre **Tracer** y carga dicho log. Comprueba que el parámetro **ESS (Effective Sample Size)** sea superior a 200 en los apartados principales para garantizar que la simulación ha convergido correctamente. Si dispones de Tracer configurado en tu entorno Linux local, simplemente ejecútalo desde tu terminal:

```bash
tracer analisis.log
```

4. **Construcción del árbol consenso (TreeAnnotator)**: Pasa el archivo `analisis.trees` por **TreeAnnotator** para resumir la colección de árboles en un único árbol de Máxima Clado de Credibilidad (MCC). Puedes ejecutarlo por línea de comandos indicando el porcentaje de descartes (*burn-in*):

```bash
treeannotator -burnin 10 analisis.trees mi_arbol_mcc.tree
```

5. **Visualización (FigTree)**: Abre `mi_arbol_mcc.tree` en **FigTree** para visualizar el árbol final con sus respectivas fechas de divergencia y soportes de nodos.

---