# Manual del Clúster HPC - Grupo de Virus Respiratorios (VirResp)

Bienvenido al repositorio central de herramientas bioinformáticas del grupo VirResp. Este sistema está diseñado para que puedas ejecutar análisis de desduplicación, alineamiento y filogenia **sin necesidad de programar** ni configurar manualmente el gestor de colas del clúster.

Tú eliges unas pocas opciones en un menú; los scripts se encargan del resto (reservar memoria, preparar el análisis, guardar los resultados).

---

## Estructura del Repositorio

El proyecto `VirResp-cluster` está organizado en dos niveles para evitar modificaciones accidentales del código de cómputo:

```
VirResp-cluster/
├── environment.yml       <-- Uso exclusivo del administrador (Germán)
├── README.md             <-- Este manual de instrucciones
├── plantillas/           <-- Scripts INTERACTIVOS que ejecuta el usuario
│   ├── lanzar_cdhit.sh
│   ├── lanzar_mafft.sh
│   ├── lanzar_fasttree.sh
│   ├── lanzar_iqtree.sh
│   └── lanzar_beast.sh
└── scripts/              <-- Scripts MAESTROS gestionados por el sistema (no tocar)
    ├── master_cdhit.sh
    ├── master_mafft.sh
    ├── master_fasttree.sh
    ├── master_iqtree.sh
    └── master_beast.sh
```

---

## Índice de Navegación Rápida

- [Mini-glosario: qué es cada cosa](#mini-glosario-qué-es-cada-cosa)
- [Guía Rápida de Trabajo (Paso a Paso)](#guía-rápida-de-trabajo-paso-a-paso)
  * [Paso 1: Localizar archivos en tu ordenador](#paso-1-localizar-archivos-en-tu-ordenador)
  * [Paso 2: Subir secuencias con rsync](#paso-2-subir-secuencias-con-rsync)
  * [Paso 3: Conectarse por SSH al clúster](#paso-3-conectarse-por-ssh-al-clúster)
  * [Paso 4: Entrar en tu directorio de trabajo](#paso-4-entrar-en-tu-directorio-de-trabajo)
  * [Paso 5: Ejecución del flujo bioinformático](#paso-5-ejecución-del-flujo-bioinformático)
  * [Paso 6: Descargar los resultados](#paso-6-descargar-los-resultados)
- [Guía Detallada de las Plantillas Interactivas](#guía-detallada-de-las-plantillas-interactivas)
  * [1. CD-HIT](#1-lanzador-de-cd-hit-lanzar_cdhitsh)
  * [2. MAFFT](#2-lanzador-de-mafft-lanzar_mafftsh)
  * [3. FastTree (árbol rápido)](#3-lanzador-de-fasttree-lanzar_fasttreesh)
  * [4. IQ-TREE 2 (árbol riguroso)](#4-lanzador-de-iq-tree-2-lanzar_iqtreesh)
  * [5. BEAST + BEAGLE](#5-lanzador-de-beast--beagle-lanzar_beastsh)
- [¿FastTree o IQ-TREE? Cuál elegir](#fasttree-o-iq-tree-cuál-elegir)
- [Qué hacer si algo falla](#qué-hacer-si-algo-falla)
- [Monitoreo y Gestión de Trabajos en Slurm](#monitoreo-y-gestión-de-trabajos-en-slurm)
- [Comandos Básicos de Linux](#comandos-básicos-de-linux)
- [Advertencia sobre la Eliminación de Archivos](#advertencia-sobre-la-eliminación-de-archivos)
- [Para el administrador](#para-el-administrador)

---

## Mini-glosario: qué es cada cosa

| Palabra | Qué significa, en cristiano |
| --- | --- |
| **Clúster** | Un conjunto de ordenadores potentes compartidos. Tú les mandas el trabajo y lo hacen ellos. |
| **Slurm / cola** | El "organizador de turnos" del clúster. Tu trabajo espera su turno y se ejecuta cuando hay hueco. |
| **FASTA** | El formato de archivo de secuencias (cada secuencia empieza por una línea con `>` y su nombre). |
| **Alineamiento** | Tus secuencias colocadas unas sobre otras para que se puedan comparar posición a posición. Es obligatorio antes de hacer un árbol. |
| **Desduplicar** | Quitar secuencias idénticas (dejando una representante) para no repetir trabajo. |
| **Árbol filogenético** | El dibujo de parentesco entre tus secuencias. Se guarda en un archivo `.treefile`. |
| **Soporte de ramas** | Un número de 0 a 1 en cada rama que indica lo fiable que es esa agrupación. |
| **Bootstrap** | Un método clásico para calcular ese soporte repitiendo el análisis muchas veces. |
| **Modelo evolutivo** | Cómo se supone que cambian las secuencias con el tiempo. Los lanzadores eligen uno adecuado por ti. |

---

## Guía Rápida de Trabajo (Paso a Paso)

El flujo completo se compone de 3 fases principales (más una fase filodinámica opcional con BEAST). Cada herramienta toma la salida del paso anterior.

```
[FASTA crudo] ──► CD-HIT ──► [FASTA desduplicado] ──► MAFFT ──► [Alineamiento FASTA] ──┬─► FastTree ──► [Árbol rápido]
                                                                                        ├─► IQ-TREE  ──► [Árbol riguroso]
                                                                                        └─► BEAST (XML) ──► [Análisis Bayesiano]
```

> **Resumen:** primero quitas repetidas (CD-HIT), luego alineas (MAFFT), y con el alineamiento haces el árbol (FastTree o IQ-TREE).

---

### Paso 1: Localizar archivos en tu ordenador

Abre la terminal en tu equipo local y navega hasta la carpeta donde tienes tu archivo de secuencias FASTA (por ejemplo, `mis_secuencias.fasta`).
> **ENTORNO: TU ORDENADOR (LOCAL)**

```
cd /ruta/a/tu/carpeta/local
ls -l mis_secuencias.fasta
```

---

### Paso 2: Subir secuencias con rsync

Envía el archivo desde tu ordenador al directorio de tu usuario dentro del clúster.
> **ENTORNO: TU ORDENADOR (LOCAL)**

```
# Sustituye 'tu_usuario' por tu nombre de usuario del clúster
rsync -avz -e "ssh -p < numero de puerto >" ./mis_secuencias.fasta tu_usuario@host:/ruta/carpeta/grupo/tu_usuario/
```

---

### Paso 3: Conectarse por SSH al clúster

Abre la conexión remota con el servidor de cómputo introduciendo tu clave personal.
> **ENTORNO: CLÚSTER (PORTUTATIS)**

```
ssh -p < numero de puerto > tu_usuario@host
```

---

### Paso 4: Entrar en tu directorio de trabajo

Dirígete a tu carpeta de usuario y verifica la presencia del archivo subido.
> **ENTORNO: CLÚSTER (PORTUTATIS)**

```
cd /ruta/carpeta/grupo/tu_usuario/
ls -l
```

---

### Paso 5: Ejecución del flujo bioinformático

Ejecuta siempre los lanzadores **desde la carpeta donde está tu archivo**: así los registros del trabajo quedan junto a tus datos.

#### 5.1. Desduplicación con CD-HIT

Elimina secuencias idénticas o redundantes antes de realizar el alineamiento.

- **Entrada:** `mis_secuencias.fasta`
- **Salida esperada:** Fichero desduplicado `mis_secuencias_representatives.fasta` e informe de haplotipos `mis_secuencias.haplotypes.txt`.

```
/ruta/carpeta/grupo/VirResp-cluster/plantillas/lanzar_cdhit.sh mis_secuencias.fasta
```

#### 5.2. Alineamiento con MAFFT

Alinea las secuencias representativas obtenidas en el paso anterior.

- **Entrada:** `mis_secuencias_representatives.fasta`
- **Salida esperada:** Fichero alineado `mis_secuencias_representatives_aligned.fasta`.

```
/ruta/carpeta/grupo/VirResp-cluster/plantillas/lanzar_mafft.sh mis_secuencias_representatives.fasta
```

#### Comprobar el estado a mitad del proceso

Para consultar si tu trabajo ha terminado o sigue en ejecución:

```
squeue -u $USER
```

Si **no aparece nada**, ha terminado (bien o mal): mira los archivos de resultados y los `.out`/`.err`.

#### 5.3. Árbol filogenético rápido con FastTree

Genera un árbol de Máxima Verosimilitud en minutos u horas, incluso con miles de secuencias. Ideal para explorar tus datos o cuando hay muchísimas secuencias.

- **Entrada:** `mis_secuencias_representatives_aligned.fasta`
- **Salida esperada:** Carpeta `mis_secuencias_representatives_aligned_fasttree_results/` con el árbol (`.treefile`), un archivo con los parámetros usados y un registro detallado.

```
/ruta/carpeta/grupo/VirResp-cluster/plantillas/lanzar_fasttree.sh mis_secuencias_representatives_aligned.fasta
```

#### 5.4. Árbol filogenético riguroso con IQ-TREE 2

Genera un árbol de Máxima Verosimilitud eligiendo el mejor modelo evolutivo. Es más lento, pero es el que se suele usar para el árbol definitivo de un artículo.

- **Entrada:** `mis_secuencias_representatives_aligned.fasta`
- **Salida esperada:** Carpeta `mis_secuencias_representatives_aligned_iqtree_results/` con el archivo del árbol `.treefile`.

```
/ruta/carpeta/grupo/VirResp-cluster/plantillas/lanzar_iqtree.sh mis_secuencias_representatives_aligned.fasta
```

#### 5.5. Análisis Filodinámico con BEAST (Opcional)

Si has preparado previamente un archivo `.xml` en BEAUti en tu ordenador local:

- **Entrada:** `mi_analisis.xml`
- **Salida esperada:** Carpeta `mi_analisis_beast_results/` con archivos `.log` y `.trees`.

> **Importante:** el archivo `.xml` generado por BEAUti ya contiene **toda la información** (secuencias alineadas + configuración del modelo). No necesitas subir el FASTA ni el NEXUS al clúster. Solo el XML.

```
/ruta/carpeta/grupo/VirResp-cluster/plantillas/lanzar_beast.sh mi_analisis.xml
```

---

### Paso 6: Descargar los resultados

Una vez finalizados los trabajos en el clúster, vuelve a la terminal de tu ordenador personal (LOCAL) para recuperar los archivos resultantes.
> **ENTORNO: TU ORDENADOR (LOCAL)**

```
# Opción A: bajar solo la carpeta de resultados de un análisis (recomendado)
rsync -avz -e "ssh -p < numero de puerto >" tu_usuario@host:/ruta/carpeta/grupo/tu_usuario/mis_secuencias_representatives_aligned_fasttree_results/ ./resultados_virresp/

# Opción B: bajar toda tu carpeta de usuario (puede ser mucho volumen)
rsync -avz -e "ssh -p < numero de puerto >" tu_usuario@host:/ruta/carpeta/grupo/tu_usuario/ ./resultados_virresp/
```

Para **ver el árbol**, abre el archivo `.treefile` con un programa como [FigTree](https://github.com/rambaut/figtree) o con la web [iTOL](https://itol.embl.de). El árbol se guarda **sin raíz**: tendrás que enraizarlo tú (con un grupo externo o por el punto medio) en el programa.

---

## Guía Detallada de las Plantillas Interactivas

Los scripts de la carpeta `plantillas/` no hacen el cálculo directamente en el nodo de entrada. En su lugar: revisan tu archivo, te muestran un menú sencillo, calculan los recursos necesarios (núcleos y memoria) y envían el trabajo a la cola (Slurm), ejecutando el script correspondiente de `scripts/master_*.sh`.

---

### 1. Lanzador de CD-HIT (`lanzar_cdhit.sh`)

- **Uso:** `/ruta/carpeta/grupo/VirResp-cluster/plantillas/lanzar_cdhit.sh <archivo.fasta>`
- **Valida e inspecciona:** Comprueba la existencia del archivo, valida que contenga cabeceras FASTA (`>`), cuenta el número total de secuencias y calcula la longitud media en posiciones.

#### Opciones interactivas que te solicitará:

1. **Tipo de Secuencia (NT vs AA):**
  - `1) Nucleótidos / ADN / ARN [Genomas, genes]`: Ejecuta el algoritmo **`cd-hit-est`**. *(El tamaño de palabra `-n` se asigna automáticamente entre 5 y 10 según el umbral).*
  - `2) Aminoácidos / Proteínas`: Ejecuta el algoritmo **`cd-hit`**. *(El tamaño de palabra `-n` se asigna automáticamente entre 2 y 5 según el umbral).*
2. **Umbral de Similitud (`-c`):**
  - `1) 100% Identidad (-c 1.0)`: Desduplicación pura (haplotipado estricto).
  - `2) 99% Identidad (-c 0.99)`: Colapsa secuencias con pequeños errores de lectura o microvariación.
  - `3) 95% Identidad (-c 0.95)`: Agrupa por linajes o subvariantes genómicas.
  - `4) Personalizado`: Permite introducir un umbral manual entre 0.80 y 1.0.
3. **Cobertura mínima de la secuencia más corta (`-aS`):** Por defecto `1.0` (100%).
4. **Resumen y confirmación:** Muestra el ejecutable seleccionado (`cd-hit` o `cd-hit-est`), los parámetros técnicos y los recursos que asignará antes de pedir confirmación (`s/n`).

#### Gestión de Recursos:

- **Datasets estándar (≤ 50.000 secuencias):** Partición `short_idx`, tiempo máx: 4 horas, 8 CPUs, 32 GB RAM.
- **Datasets masivos (> 50.000 secuencias):** Partición `middle_idx`, tiempo máx: 24 horas, 16 CPUs, 64 GB RAM.

---

### 2. Lanzador de MAFFT (`lanzar_mafft.sh`)

- **Uso:** `/ruta/carpeta/grupo/VirResp-cluster/plantillas/lanzar_mafft.sh <archivo.fasta>`
- **Valida e inspecciona:** Comprueba el formato FASTA, cuenta el número de muestras y evalúa la longitud media.

#### Opciones interactivas que te solicitará:

1. **Estrategias de Alineamiento:**
  - `1) AUTO`: Recomendado por defecto. Selecciona automáticamente el mejor algoritmo.
  - `2) L-INS-i`: Alta precisión (`--maxiterate 1000 --localpair`). Ideal para < 200 secuencias o regiones muy variables.
  - `3) FFT-NS-2`: Método rápido de 2 etapas (`--retree 2`). Adecuado para miles de secuencias.
  - `4) PartTree`: Algoritmo ultra-masivo (`--parttree --retree 1`). Obligatorio si tienes > 10.000 genomas completos.
2. **Ajuste de orientación de hebras (`--adjustdirection`):** Activo por defecto (`S/n`) para detectar secuencias en sentido inverso.
3. **Reordenar secuencias por similitud (`--reorder`):** Inactivo por defecto (`s/N`).
4. **Penalización de apertura de GAPs (`--op`):** Valor predeterminado `1.53`.

#### Gestión de Recursos:

- **Estándar / AUTO:** Partición `short_idx`, 12 horas, 16 CPUs, 64 GB RAM.
- **Estrategia L-INS-i:** Partición `middle_idx`, 24 horas, 16 CPUs, 64 GB RAM.
- **Estrategia PartTree o datasets > 10.000 secuencias:** Partición `long_idx`, 5 días, 32 CPUs, 300 GB RAM.

---

### 3. Lanzador de FastTree (`lanzar_fasttree.sh`)

- **Uso:** `/ruta/carpeta/grupo/VirResp-cluster/plantillas/lanzar_fasttree.sh <alineamiento.fasta>`
- **Necesitas:** un FASTA **ya alineado** (salida de MAFFT).

#### Qué comprueba antes de preguntarte nada

El script revisa tu archivo y **te avisa en lenguaje claro** si hay un problema, sin llegar a gastar tiempo de cálculo:

- Que sea realmente un FASTA y tenga al menos 4 secuencias.
- Que **todas las secuencias midan lo mismo** (si no, no están alineadas: te manda a MAFFT).
- Que **no haya nombres repetidos**.
- Que los nombres no tengan símbolos que estropean el árbol: `( ) , : ; [ ]`.
- Si son nucleótidos o proteínas (lo detecta él solo y te lo propone).

Tu archivo original **no se modifica**.

#### Las 3 preguntas que te hará

1. **¿Qué tipo de secuencias son?** El script te propone la opción correcta; normalmente solo tienes que pulsar Intro.
2. **¿Qué tipo de búsqueda quieres?**
  - `1) Estándar` (recomendado): rápido y fiable para casi todos los casos.
  - `2) Cuidadoso`: dedica más tiempo a buscar el mejor árbol posible (varias veces más lento). Para el árbol definitivo.
3. **¿Cómo medir la fiabilidad de las ramas?**
  - `1) Rápida (SH-like)` (recomendado): se calcula sin coste extra. **No es un bootstrap**, no la llames así en un artículo.
  - `2) Bootstrap (100 repeticiones)`: el método clásico; más lento. Solo disponible hasta 5.000 secuencias.
  - `3) Sin valores`: lo más rápido.

Después muestra un resumen y pide confirmación (`s/n`) antes de enviar nada.

#### Lo que decide solo (no tienes que preocuparte)

- **Modelo evolutivo:** GTR para nucleótidos, LG para proteínas. Son los estándar para este tipo de datos.
- **Longitudes de rama refinadas** (opción `-gamma`) siempre activadas.
- **Semilla fija** (1253): si repites el análisis con el mismo archivo y opciones, sale el mismo árbol.
- **Más de 50.000 secuencias:** activa un modo de ahorro de memoria y tiempo.

#### Qué recursos reserva

Se calculan según **nº de secuencias × longitud del alineamiento** y según lo que elijas (cuidadoso y bootstrap necesitan más):

| Tamaño del trabajo | Cola | Tiempo máx. | Núcleos |
| --- | --- | --- | --- |
| Pequeño (hasta ~5.000 secuencias y poco volumen total) | `short_idx` | 12 horas | 4 |
| Mediano | `middle_idx` | 48 horas | 8 |
| Grande (> 50.000 secuencias o muchísimo volumen) | `long_idx` | 5 días | 16 |

La memoria se estima de forma orientativa (mínimo 16 GB, máximo 300 GB). Si un trabajo se queda sin memoria o sin tiempo, consulta [Qué hacer si algo falla](#qué-hacer-si-algo-falla).

#### Qué obtienes

Una carpeta `<tu_archivo>_fasttree_results/` con tres archivos, todos con la fecha y hora del análisis en el nombre (**así nunca se sobrescriben unos análisis con otros**):

- `*.treefile`: el árbol (sin raíz). Los números sobre las ramas van de 0 a 1 (0,95 = 95 %).
- `*_parametros.txt`: todo lo necesario para reproducir y describir el análisis (parámetros, versión del programa, nº de secuencias, huella del archivo de entrada, cita).
- `*.log`: registro detallado del programa.

#### Cómo leer el soporte de ramas (orientativo)

Los valores SH-like cercanos a 1 indican agrupaciones muy fiables; los valores bajos indican que esa rama es dudosa. Como orientación, muchos grupos consideran fiables las ramas por encima de ~0,9, pero no son equivalentes a un bootstrap: si necesitas valores publicables, usa la opción de bootstrap o IQ-TREE.

#### Cómo citar

Price MN, Dehal PS, Arkin AP (2010) *FastTree 2 – Approximately Maximum-Likelihood Trees for Large Alignments.* PLoS ONE 5(3): e9490.

---

### 4. Lanzador de IQ-TREE 2 (`lanzar_iqtree.sh`)

- **Uso:** `/ruta/carpeta/grupo/VirResp-cluster/plantillas/lanzar_iqtree.sh <alineamiento.fasta>`
- **Valida e inspecciona:** Comprueba la validez del alineamiento y calcula la longitud exacta del alineamiento incluyendo gaps/N.

#### Opciones interactivas que te solicitará:

1. **Modelo de Sustitución Evolutiva:**
  - `1) Búsqueda automática con ModelFinder (-m MF)` [Recomendado].
  - `2) GTR+F+I+G4`: Modelo nucleotídico completo y robusto.
  - `3) GTR+G`: Modelo estándar de Máxima Verosimilitud.
  - `4) HKY85+G`: Para secuencias poco divergentes.
  - `5) TN93+G`: Modelo de Tamura-Nei.
  - `6) LG+F+G4`: Modelo estándar para alineamientos de aminoácidos / proteínas.
  - `7) Personalizado`: Escribir manualmente la sintaxis exacta de IQ-TREE (ej. `GTR+F+R4`).
2. **Método de Soporte de Ramas (Bootstrap):**
  - `1) Ultrafast Bootstrap (UFBoot)`: 1000 réplicas (`-bb 1000`).
  - `2) SH-aLRT test + UFBoot`: Doble validación con 1000 réplicas cada una (`-alrt 1000 -bb 1000`).
  - `3) Sin Bootstrap`: Genera solo la topología del árbol.

#### Gestión de Recursos:

- **Alineamientos pequeños/Modelos fijos:** Partición `short_idx`, 12 horas, 16 CPUs, 64 GB RAM.
- **Con ModelFinder (MF) o > 1.000 secuencias:** Partición `middle_idx`, 48 horas, 32 CPUs, 128 GB RAM.
- **Datasets masivos (> 5.000 secuencias):** Escalado automático a partición `long_idx`, 5 días, 32 CPUs, 300 GB RAM (`-mem 300G`).

---

### 5. Lanzador de BEAST + BEAGLE (`lanzar_beast.sh`)

- **Uso:** `/ruta/carpeta/grupo/VirResp-cluster/plantillas/lanzar_beast.sh <analisis.xml>`
- **Valida e inspecciona automáticamente:** Comprueba que el XML existe, elimina saltos de línea incompatibles de Windows (`\r`), extrae la longitud de la cadena MCMC (`chainLength`) y calcula el número de taxones/secuencias.

#### Configuración y asignación automática de recursos:

Este script **es 100% automático** y no requiere seleccionar opciones manuales. Aplica la configuración más estable para evitar errores de desbordamiento numérico:

- **Motor de cálculo:** CPU Multinúcleo + BEAGLE SSE con Doble Precisión (`-beagle -beagle_CPU -beagle_SSE -beagle_double`).
- **Asignación inteligente:**
  * **Análisis Estándar (`chainLength` ≤ 50M y ≤ 300 taxones):** Partición `middle_idx`, tiempo máx: 48 horas, 16 CPUs, 64 GB RAM.
  * **Datasets Masivos o MCMC Larga (`chainLength` > 50M o > 300 taxones):** Escalado automático a partición `long_idx`, tiempo máx: 5 días, 32 CPUs, 128 GB RAM.

---

## FastTree o IQ-TREE: cuál elegir

Ambos hacen el mismo tipo de árbol (máxima verosimilitud). La diferencia es **velocidad frente a rigor**.

| | FastTree | IQ-TREE |
| --- | --- | --- |
| **Velocidad** | Muy rápido | Más lento |
| **Muchas secuencias (miles)** | Muy adecuado | Puede tardar mucho |
| **Modelo evolutivo** | Fijo (GTR o LG), sin comparar modelos | Busca automáticamente el mejor (ModelFinder) |
| **Fiabilidad de ramas** | SH-like (rápida) o bootstrap lento | UFBoot y SH-aLRT (estándar en publicaciones) |
| **Úsalo para** | Explorar, comprobar que el alineamiento tiene sentido, conjuntos muy grandes | El árbol definitivo de un artículo |

**Recomendación general (orientativa):** empieza con FastTree para echar un primer vistazo a tus datos; cuando tengas el conjunto definitivo, haz el árbol final con IQ-TREE. Si tienes decenas de miles de secuencias, FastTree suele ser la única opción práctica.

---

## Qué hacer si algo falla

Primero mira el final de los archivos de registro de tu trabajo (están en la misma carpeta que tu archivo):

```
tail -n 30 fasttree_12345.err
tail -n 30 fasttree_12345.out
```

> El archivo `.err` de FastTree **siempre tiene mensajes de progreso**. Que tenga texto no significa que haya un error.

| Qué ves | Qué significa | Qué hacer |
| --- | --- | --- |
| Antes de enviar: *"Tus secuencias NO miden todas lo mismo"* | El archivo no está alineado | Alinéalo primero con `lanzar_mafft.sh` |
| Antes de enviar: *"secuencias con el MISMO nombre"* | Hay nombres de secuencia repetidos | Renombra las repetidas en tu FASTA |
| Aviso de símbolos `( ) , : ; [ ]` en los nombres | Pueden estropear el árbol | Cámbialos por guiones bajos y relanza |
| `squeue` muestra `PD` mucho tiempo | Esperando turno en la cola | Es normal; si pasan días, avisa al administrador |
| El trabajo desaparece sin resultados y en el `.err` pone `oom-kill` o `out of memory` | Se quedó sin memoria | Avisa al administrador indicando el número de trabajo |
| El trabajo desaparece sin resultados y en el `.err` pone `TIME LIMIT` o `DUE TO TIME LIMIT` | Se quedó sin tiempo | Avisa al administrador indicando el número de trabajo |
| Mensaje `ERROR CRÍTICO` en el `.out` o `.err` | Falló FastTree o el entorno | Copia las últimas 30 líneas del `.err` y envíaselas al administrador |

Para consultar cuánto usó un trabajo ya terminado (útil para el administrador):

```
sacct -j 12345 --format=JobID,State,Elapsed,MaxRSS,ReqMem
```

Si un análisis falla a mitad, el script guarda el registro parcial en la carpeta de resultados con la palabra `FALLO` en el nombre, para ayudar a diagnosticar.

---

## Monitoreo y Gestión de Trabajos en Slurm

Los lanzadores interactivos guardan los archivos de salida y de error en el mismo directorio donde está tu archivo de entrada (el FASTA o el XML que le pasas como argumento). Por eso se recomienda ejecutar el script desde esa misma carpeta: así los registros quedan junto a tus datos y resultados.

### Archivos de salida de Slurm generados:

- `cdhit_<JOB_ID>.out` / `cdhit_<JOB_ID>.err`
- `mafft_<JOB_ID>.out` / `mafft_<JOB_ID>.err`
- `fasttree_<JOB_ID>.out` / `fasttree_<JOB_ID>.err`
- `iqtree_<JOB_ID>.out` / `iqtree_<JOB_ID>.err`
- `beast_<JOB_ID>.out` / `beast_<JOB_ID>.err`

(`<JOB_ID>` es el número de trabajo que te muestra el lanzador al enviarlo.)

### Comandos de control en la terminal:

```
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

- **`R` (Running):** El trabajo se está ejecutando activamente en los nodos de cálculo.
- **`PD` (Pending):** El trabajo está en cola esperando que queden libres los núcleos/memoria solicitados.
- **`CG` (Completing):** El trabajo ha terminado y está guardando los archivos finales en el disco.

---

## Comandos Básicos de Linux

| Comando                    | Función / Uso                                                  |
| -------------------------- | -------------------------------------------------------------- |
| `cd /ruta/`                | Cambiar de carpeta/directorio.                                 |
| `cd ..`                    | Volver a la carpeta anterior (subir un nivel).                 |
| `ls -lh`                   | Listar archivos mostrando el tamaño de forma legible (MB, GB). |
| `pwd`                      | Mostrar la ruta completa del directorio donde te encuentras.   |
| `mkdir nueva_carpeta`      | Crear una carpeta nueva.                                       |
| `cp origen destino`        | Copiar un archivo a otra ubicación.                            |
| `mv origen destino`        | Mover o cambiar de nombre un archivo o carpeta.                |
| `head -n 20 archivo.fasta` | Ver las primeras 20 líneas de un archivo.                      |
| `tail -n 20 archivo.out`   | Ver las últimas 20 líneas de un archivo de salida o log.       |
| `grep -c ">" archivo.fasta`| Contar cuántas secuencias tiene un FASTA.                      |

---

## Advertencia sobre la Eliminación de Archivos

> **¡ATENCIÓN! EN EL SISTEMA DE ARCHIVOS DEL CLÚSTER NO EXISTE PAPELERA DE RECICLAJE**
>
> - **`rm archivo`**: Borra un archivo de forma definitiva.
> - **`rm -r carpeta`**: Borra una carpeta completa y todo su contenido de forma irreversible.
>
> Revisa cuidadosamente la ruta y el nombre del archivo antes de presionar `Enter` al ejecutar comandos de borrado.

---

## Para el administrador

Esta sección es solo para quien mantiene el repositorio.

**Entorno conda.** Se define en `environment.yml` (instrucciones de creación, actualización y congelado en la cabecera del archivo). Tras cada cambio, genera y sube `environment.lock.yml` con las versiones exactas.

**Configuración opcional de `master_fasttree.sh`** (variables de entorno):

- `VIRRESP_ENV`: ruta alternativa del entorno conda.
- `FASTTREE_BIN`: fuerza un ejecutable concreto. Por ejemplo, `FASTTREE_BIN=FastTreeDbl ./lanzar_fasttree.sh x.fasta` usa doble precisión, recomendable si hay conjuntos con secuencias casi idénticas sin desduplicar (esta versión usa un solo núcleo).

**Cómo decide FastTree los recursos.** `lanzar_fasttree.sh` estima el trabajo como secuencias × longitud (multiplicado por 3 en modo cuidadoso y por 4 con bootstrap), y la memoria como 3 × secuencias × longitud × 16 B (nucleótidos) u 80 B (proteínas), +5 GB, redondeada a múltiplos de 8 GB (16-300 GB). Son estimaciones iniciales: contrástalas con `sacct ... MaxRSS` en los primeros trabajos reales y ajusta los umbrales.

**Revisión del código.** Antes de subir cambios:

```
shellcheck plantillas/*.sh scripts/*.sh
```

**Citas de las herramientas** (para métodos de los artículos): FastTree 2 (Price et al., 2010, PLoS ONE 5:e9490); consulta `*_parametros.txt` de cada análisis para versiones exactas.
