# Manual del Clúster HPC - Grupo de Virus Respiratorios (VirResp)

Bienvenido al repositorio central de herramientas bioinformáticas del grupo VirResp. Este sistema está diseñado para que puedas ejecutar análisis de alineamiento, desduplicación y filogenia sin necesidad de programar ni configurar manualmente el gestor de colas del clúster.

---

## Estructura del Repositorio

El proyecto `VirResp-cluster` está organizado en dos niveles para evitar modificaciones accidental del código de cómputo:


```text
VirResp-cluster/
├── environment.yml       <-- Uso exclusivo del administrador (Germán)
├── README.md             <-- Este manual de instrucciones
├── plantillas/           <-- Scripts INTERACTIVOS que ejecuta el usuario
│   ├── lanzar_mafft.sh
│   ├── lanzar_cdhit.sh
│   ├── lanzar_iqtree.sh
│   └── lanzar_beast.sh
└── scripts/              <-- Scripts MAESTROS gestionados por el sistema
    ├── master_mafft.sh
    ├── master_cdhit.sh
    ├── master_iqtree.sh
    └── master_beast.sh
```


---

## Índice de Navegación Rápida

* [🎯 Guía Rápida de Trabajo (Paso a Paso)](#-guía-rápida-de-trabajo-paso-a-paso)
  * [Paso 1: Localizar archivos en tu ordenador](#paso-1-localizar-archivos-en-tu-ordenador)
  * [Paso 2: Subir secuencias con rsync](#paso-2-subir-secuencias-con-rsync)
  * [Paso 3: Conectarse por SSH al clúster](#paso-3-conectarse-por-ssh-al-clúster)
  * [Paso 4: Entrar en tu directorio de trabajo](#paso-4-entrar-en-tu-directorio-de-trabajo)
  * [Paso 5: Ejecución del flujo bioinformático](#paso-5-ejecución-del-flujo-bioinformático)
  * [Paso 6: Descargar los resultados](#paso-6-descargar-los-resultados)
* [🛠️ Guía Detallada de las Plantillas Interactivas](#️-guía-detallada-de-las-plantillas-interactivas)
  * [1. Lanzador de CD-HIT (`lanzar_cdhit.sh`)](#1-lanzador-de-cd-hit-lanzar_cdhitsh)
  * [2. Lanzador de MAFFT (`lanzar_mafft.sh`)](#2-lanzador-de-mafft-lanzar_mafftsh)
  * [3. Lanzador de IQ-TREE 2 (`lanzar_iqtree.sh`)](#3-lanzador-de-iq-tree-2-lanzar_iqtreesh)
  * [4. Lanzador de BEAST + BEAGLE (`lanzar_beast.sh`)](#4-lanzador-de-beast--beagle-lanzar_beastsh)
* [📊 Monitoreo y Gestión de Trabajos en Slurm](#-monitoreo-y-gestión-de-trabajos-en-slurm)
* [🛠️ Comandos Básicos de Linux](#️-comandos-básicos-de-linux)
* [⚠️ Advertencia sobre la Eliminación de Archivos](#️-advertencia-sobre-la-eliminación-de-archivos)

---

## 🎯 Guía Rápida de Trabajo (Paso a Paso)

El flujo de trabajo completo se compone de 3 fases principales (más una fase filodinámica opcional con BEAST). Cada herramienta toma la salida del paso anterior.


```text
[FASTA crudo] ──► CD-HIT ──► [FASTA desduplicado] ──► MAFFT ──► [Alineamiento FASTA] ──► IQ-TREE ──► [Árbol ML]
                                                                                           └──► BEAST (XML) ──► [Análisis Bayesiano]
```


---

### Paso 1: Localizar archivos en tu ordenador
Abre la terminal en tu equipo local y navega hasta la carpeta donde tienes tu archivo de secuencias FASTA (por ejemplo, `mis_secuencias.fasta`).

> 🟦 **ENTORNO: TU ORDENADOR (LOCAL)**


```bash
cd /ruta/a/tu/carpeta/local
ls -l mis_secuencias.fasta
```


---

### Paso 2: Subir secuencias con rsync
Envía el archivo desde tu ordenador al directorio de tu usuario dentro del clúster.

> 🟦 **ENTORNO: TU ORDENADOR (LOCAL)**


```bash
# Sustituye 'tu_usuario' por tu nombre de usuario del clúster
rsync -avz -e "ssh -p < numero de puerto >" ./mis_secuencias.fasta tu_usuario@host:/ruta/carpeta/grupo/tu_usuario/
```


---

### Paso 3: Conectarse por SSH al clúster
Abre la conexión remota con el servidor de cómputo introduciendo tu clave personal.

> 🟩 **ENTORNO: CLÚSTER (PORTUTATIS)**


```bash
ssh -p < numero de puerto > tu_usuario@host
```


---

### Paso 4: Entrar en tu directorio de trabajo
Dirígete a tu carpeta de usuario y verifica la presencia del archivo subido.

> 🟩 **ENTORNO: CLÚSTER (PORTUTATIS)**


```bash
cd /ruta/carpeta/grupo/tu_usuario/
ls -l
```


---

### Paso 5: Ejecución del flujo bioinformático

#### 5.1. Desduplicación con CD-HIT
Elimina secuencias idénticas o redundantes antes de realizar el alineamiento.

* **Entrada:** `mis_secuencias.fasta`
* **Salida esperada:** Fichero desduplicado `mis_secuencias_representatives.fasta` e informe de haplotipos `mis_secuencias.haplotypes.txt`.


```bash
/ruta/carpeta/grupo/VirResp-cluster/plantillas/lanzar_cdhit.sh mis_secuencias.fasta
```


#### 5.2. Alineamiento con MAFFT
Alinea las secuencias representativas obtenidas en el paso anterior.

* **Entrada:** `mis_secuencias_representatives.fasta`
* **Salida esperada:** Fichero alineado `mis_secuencias_representatives_aligned.fasta`.


```bash
/ruta/carpeta/grupo/VirResp-cluster/plantillas/lanzar_mafft.sh mis_secuencias_representatives.fasta
```


#### 💡 Comprobar el estado en Slurm a mitad del proceso
Para consultar si tu trabajo en MAFFT o CD-HIT ha terminado o sigue en ejecución:


```bash
squeue -u $USER
```


#### 5.3. Reconstrucción Filogenética con IQ-TREE 2
Genera un árbol de Máxima Verosimilitud (ML) a partir del alineamiento.

* **Entrada:** `mis_secuencias_representatives_aligned.fasta`
* **Salida esperada:** Carpeta `mis_secuencias_representatives_aligned_iqtree_results/` con el archivo del árbol `.treefile`.


```bash
/ruta/carpeta/grupo/VirResp-cluster/plantillas/lanzar_iqtree.sh mis_secuencias_representatives_aligned.fasta
```


#### 5.4. Análisis Filodinámico con BEAST (Opcional)
Si has preparado previamente un archivo `.xml` en BEAUti en tu ordenador local:

* **Entrada:** `mi_analisis.xml`
* **Salida esperada:** Carpeta `mi_analisis_beast_results/` con archivos `.log` y `.trees`.

> ⚠️ **Importante:** el archivo `.xml` generado por BEAUti ya contiene **toda la información** (secuencias alineadas + configuración del modelo). No necesitas subir el FASTA ni el NEXUS al clúster. Solo el XML.


```bash
/ruta/carpeta/grupo/VirResp-cluster/plantillas/lanzar_beast.sh mi_analisis.xml
```

---

### Paso 6: Descargar los resultados
Una vez finalizados los trabajos en el clúster, vuelve a la terminal de tu ordenador personal (LOCAL) para recuperar los archivos resultantes.

> 🟦 **ENTORNO: TU ORDENADOR (LOCAL)**


```bash
rsync -avz -e "ssh -p < numero de puerto >" tu_usuario@host:/ruta/carpeta/grupo/tu_usuario/ ./resultados_virresp/
```


---

## 🛠️ Guía Detallada de las Plantillas Interactivas

Los scripts de la carpeta `plantillas/` no realizan el cómputo directamente en el nodo de login. En su lugar, realizan un análisis del archivo de entrada, muestran un menú interactivo en la terminal para que selecciones los parámetros deseados, calculan los recursos necesarios de CPU/RAM y envían la orden al gestor de colas Slurm (ejecutando el script correspondiente en `scripts/master_*.sh`).

---

### 1. Lanzador de CD-HIT (`lanzar_cdhit.sh`)

* **Uso:** `/ruta/carpeta/grupo/VirResp-cluster/plantillas/lanzar_cdhit.sh <archivo.fasta>`
* **Valida e inspecciona:** Comprueba la existencia del archivo, valida que contenga cabeceras FASTA (`>`), cuenta el número total de secuencias y calcula la longitud media en pares de bases (bp).

#### Opciones interactivas que te solicitará:
1. **Umbral de Similitud (`-c`):**
   * `1) 100% Identidad (-c 1.0)`: Desduplicación pura (haplotipado estricto).
   * `2) 99% Identidad (-c 0.99)`: Colapsa secuencias con pequeños errores de lectura o microvariación.
   * `3) 95% Identidad (-c 0.95)`: Agrupa por linajes o subvariantes genómicas.
   * `4) Personalizado`: Permite introducir un umbral manual entre 0.80 y 1.0.
   *(Nota: El tamaño de palabra `-n` se asigna automáticamente según el umbral: `-n 5` para <0.88, `-n 6` para <0.90, `-n 8` para <0.92 y `-n 10` para ≥0.92).*
2. **Cobertura mínima de la secuencia más corta (`-aS`):** Por defecto `1.0` (100%).
3. **Resumen y confirmación:** Te muestra los recursos Slurm que asignará y pide confirmación (`s/n`) para enviar.

#### Gestión de Recursos en Slurm:
* **Datasets estándar (≤ 50.000 secuencias):** Partición `short_idx`, tiempo máx: 4 horas, 8 CPUs, 32 GB RAM.
* **Datasets masivos (> 50.000 secuencias):** Partición `middle_idx`, tiempo máx: 24 horas, 16 CPUs, 64 GB RAM.

---

### 2. Lanzador de MAFFT (`lanzar_mafft.sh`)

* **Uso:** `/ruta/carpeta/grupo/VirResp-cluster/plantillas/lanzar_mafft.sh <archivo.fasta>`
* **Valida e inspecciona:** Comprueba el formato FASTA, cuenta el número de muestras y evalúa la longitud media.

#### Opciones interactivas que te solicitará:
1. **Estrategias de Alineamiento:**
   * `1) AUTO`: Recomendado por defecto. Selecciona automáticamente el mejor algoritmo.
   * `2) L-INS-i`: Alta precisión (`--maxiterate 1000 --localpair`). Ideal para < 200 secuencias o regiones muy variables.
   * `3) FFT-NS-2`: Método rápido de 2 etapas (`--retree 2`). Adecuado para miles de secuencias.
   * `4) PartTree`: Algoritmo ultra-masivo (`--parttree --retree 1`). Obligatorio si tienes > 10.000 genomas completos.
2. **Ajuste de orientación de hebras (`--adjustdirection`):** Activo por defecto (`S/n`) para detectar secuencias en sentido inverso.
3. **Reordenar secuencias por similitud (`--reorder`):** Inactivo por defecto (`s/N`).
4. **Penalización de apertura de GAPs (`--op`):** Valor predeterminado `1.53`.

#### Gestión de Recursos en Slurm:
* **Estándar / AUTO:** Partición `short_idx`, 12 horas, 16 CPUs, 64 GB RAM.
* **Estrategia L-INS-i:** Partición `middle_idx`, 24 horas, 16 CPUs, 64 GB RAM.
* **Estrategia PartTree o datasets > 10.000 secuencias:** Partición `long_idx`, 5 días, 32 CPUs, 300 GB RAM.

---

### 3. Lanzador de IQ-TREE 2 (`lanzar_iqtree.sh`)

* **Uso:** `/ruta/carpeta/grupo/VirResp-cluster/plantillas/lanzar_iqtree.sh <alineamiento.fasta>`
* **Valida e inspecciona:** Comprueba la validez del alineamiento y calcula la longitud exacta del alineamiento incluyendo gaps/N.

#### Opciones interactivas que te solicitará:
1. **Modelo de Sustitución Evolutiva:**
   * `1) Búsqueda automática con ModelFinder (-m MF)` [Recomendado].
   * `2) GTR+F+I+G4`: Modelo nucleotídico completo y robusto.
   * `3) GTR+G`: Modelo estándar de Máxima Verosimilitud.
   * `4) HKY85+G`: Para secuencias poco divergentes.
   * `5) TN93+G`: Modelo de Tamura-Nei.
   * `6) LG+F+G4`: Modelo estándar para alineamientos de aminoácidos / proteínas.
   * `7) Personalizado`: Escribir manualmente la sintaxis exacta de IQ-TREE (ej. `GTR+F+R4`).
2. **Método de Soporte de Ramas (Bootstrap):**
   * `1) Ultrafast Bootstrap (UFBoot)`: 1000 réplicas (`-bb 1000`).
   * `2) SH-aLRT test + UFBoot`: Doble validación con 1000 réplicas cada una (`-alrt 1000 -bb 1000`).
   * `3) Sin Bootstrap`: Genera solo la topología del árbol.

#### Gestión de Recursos en Slurm:
* **Alineamientos pequeños/Modelos fijos:** Partición `short_idx`, 12 horas, 16 CPUs, 64 GB RAM.
* **Con ModelFinder (MF) o > 1.000 secuencias:** Partición `middle_idx`, 48 horas, 32 CPUs, 128 GB RAM.
* **Datasets masivos (> 5.000 secuencias):** Escalado automático a partición `long_idx`, 5 días, 32 CPUs, 300 GB RAM (`-mem 300G`).

---

### 4. Lanzador de BEAST + BEAGLE (`lanzar_beast.sh`)

* **Uso:** `/ruta/carpeta/grupo/VirResp-cluster/plantillas/lanzar_beast.sh <analisis.xml>`
* **Valida e inspecciona:** Comprueba que el XML existe, elimina saltos de línea incompatibles de Windows (`\r`) y extrae automáticamente la longitud de la cadena MCMC (`chainLength`) mediante expresiones regulares.

#### Opciones interactivas que te solicitará:
1. **Motor de Aceleración Computacional:**
   * `1) CPU Multinúcleo + BEAGLE` (`-beagle -beagle_CPU -beagle_SSE`) [Recomendado].
   * `2) GPU + BEAGLE` (`-beagle -beagle_GPU`). Solicita aceleración por tarjeta gráfica en Slurm (`--gres=gpu:1`).
   * `3) BEAST Nativo`: Sin aceleración BEAGLE.
2. **Precisión Numérica de Cálculo:**
   * `1) Doble Precisión (-beagle_double)`: Estabilidad en filodinámica y relojes complejos.
   * `2) Precisión Simple (-beagle_single)`: Más rápido.
   *(Nota: Si seleccionas GPU, el script fuerza automáticamente `-beagle_single` para garantizar compatibilidad).*
3. **Estimación de tiempo de ejecución (para modo CPU):**
   * `1) MCMC Corto/Medio`: Hasta 48 horas en partición `middle_idx` (16 CPUs, 64 GB RAM).
   * `2) MCMC Largo/Muy pesado`: Hasta 5 días en partición `long_idx` (32 CPUs, 128 GB RAM).

---

## 📊 Monitoreo y Gestión de Trabajos en Slurm

Los lanzadores interactivos guardan los archivos de salida y de error en el mismo directorio donde está tu archivo de entrada (el FASTA o el XML que le pasas como argumento). Por eso se recomienda ejecutar el script desde esa misma carpeta: así los logs quedan junto a tus datos y resultados.

### Archivos de salida de Slurm generados:
* `cdhit_<JOB_ID>.out` / `cdhit_<JOB_ID>.err`
* `mafft_<JOB_ID>.out` / `mafft_<JOB_ID>.err`
* `iqtree_<JOB_ID>.out` / `iqtree_<JOB_ID>.err`
* `beast_<JOB_ID>.out` / `beast_<JOB_ID>.err`

### Comandos de control en la terminal:


```bash
# Ver el estado de todos tus trabajos en la cola
squeue -u $USER

# Monitorear la salida en tiempo real de un proceso en ejecución
tail -f mafft_12345.out

# Ver información detallada de un trabajo en específico
scontrol show job 12345

# Cancelar o eliminar un trabajo de la cola
scancel 12345
```


#### Significado de los estados en `squeue`:
* **`R` (Running):** El trabajo se está ejecutando activamente en los nodos de cálculo.
* **`PD` (Pending):** El trabajo está en cola esperando que queden libres las CPUs/RAM solicitadas.
* **`CG` (Completing):** El trabajo ha terminado y está guardando los archivos finales en el disco.

---

## 🛠️ Comandos Básicos de Linux

| Comando | Función / Uso |
| :--- | :--- |
| `cd /ruta/` | Cambiar de carpeta/directorio. |
| `cd ..` | Volver a la carpeta anterior (subir un nivel). |
| `ls -lh` | Listar archivos mostrando el tamaño de forma legible (MB, GB). |
| `pwd` | Mostrar la ruta completa del directorio donde te encuentras. |
| `mkdir nueva_carpeta` | Crear una carpeta nueva. |
| `cp origen destino` | Copiar un archivo a otra ubicación. |
| `mv origen destino` | Mover o cambiar de nombre un archivo o carpeta. |
| `head -n 20 archivo.fasta` | Ver las primeras 20 líneas de un archivo. |
| `tail -n 20 archivo.out` | Ver las últimas 20 líneas de un archivo de salida o log. |

---

## ⚠️ Advertencia sobre la Eliminación de Archivos

> 🛑 **¡ATENCIÓN! EN EL SISTEMA DE ARCHIVOS DEL CLÚSTER NO EXISTE PAPELERA DE RECICLAJE**
> 
> * **`rm archivo`**: Borra un archivo de forma definitiva.
> * **`rm -r carpeta`**: Borra una carpeta completa y todo su contenido de forma irreversible.
> 
> Revisa cuidadosamente la ruta y el nombre del archivo antes de presionar `Enter` al ejecutar comandos de borrado.