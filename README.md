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

## Mini-glosario del Clúster y Herramientas.

| Término | Definición / Contexto |
| --- | --- |
| **Clúster** | Conjunto de nodos de cómputo de alto rendimiento conectados en red. |
| **Slurm / Cola** | Gestor de recursos y trabajos que asigna tiempo de CPU, memoria y prioridad de ejecución. |
| **CD-HIT / CD-HIT-EST** | Herramienta ultrarrápida de agrupamiento de secuencias. `cd-hit` se emplea para proteínas y `cd-hit-est` para secuencias nucleotídicas (ADN/ARN). |
| **Haplotipado / Desduplicación** | Proceso de colapsar secuencias 100% idénticas en una única secuencia representativa para eliminar redundancia antes de análisis filogenéticos complejos. |
| **Umbral de Similitud (`-c`)** | Porcentaje de identidad de secuencia (ej. `1.0` para 100%, `0.95` para 95%) requerido para agrupar dos secuencias dentro del mismo clúster. |
| **Word Size (`-n`)** | Tamaño del oligómero corto (k-mero) que utiliza CD-HIT para el indexado preliminar. El script ajusta automáticamente este valor según el umbral `-c` seleccionado. |
| **ModelFinder (MF / MFP)** | Algoritmo integrado en IQ-TREE 2 que evalúa exhaustivamente cientos de modelos de sustitución sustitutivos (y tasas de heterogeneidad entre sitios como FreeRate `+R`) seleccionando el óptimo mediante AIC, AICc y BIC. |
| **Ultrafast Bootstrap (UFBoot)** | Algoritmo de remuestreo bootstrap aproximado que reduce el tiempo de cómputo en varios órdenes de magnitud con respecto al bootstrap estándar no paramétrico. Valores $\ge 95\%$ indican soporte fuerte. |
| **SH-aLRT** | Prueba de verosimilitud local Shimodaira-Hasegawa aproximada. Evalúa la hipótesis nula de que la longitud de la rama interna es cero. Valores $\ge 80\%$ indican ramas significativamente respaldadas. |
| **Optimización BNNI** | Optimización de intercambios de vecinos más cercanos (NNI) en los árboles resueltos por UFBoot para evitar la sobreestimación de valores de soporte bajo violaciones severas del modelo. |
| **Modelo CAT (-cat)** | Aproximación rápida a la heterogeneidad de tasas entre sitios mediante $N$ categorías discretas de velocidad de sustitución. |
| **Reoptimización Gamma (-gamma)** | Recálculo de las longitudes de rama y verosimilitud (lnL) bajo una distribución Gamma tras resolver la topología con CAT. |
| **L-INS-i / G-INS-i** | Algoritmos de alineamiento iterativo en MAFFT. `L-INS-i` optimiza la homología local (ideal para regiones variables) y `G-INS-i` la homología global a lo largo de toda la secuencia. |
| **Leave gappy regions** | Algoritmo de MAFFT que evita forzar el alineamiento en regiones no homólogas, inserciones virales únicas o extremos incompletos, dejándolos en bloques no alineados (gaps). |
| **Unalignlevel** | Parámetro numérico (0.0 a 0.8) que define el umbral de tolerancia para declarar una región como no homóloga antes de aplicar `Leave gappy regions`. |
| **BEAST** | Plataforma de análisis filogenético bayesiano mediante cadenas MCMC para la estimación de relojes moleculares, tasas de sustitución y parámetros filodinámicos de poblaciones. |
| **BEAGLE** | Librería de cómputo de alto rendimiento integrada en BEAST para la paralelización y aceleración del cálculo de verosimilitud en vectores CPU/GPU. |
| **MCMC / chainLength** | *Markov Chain Monte Carlo*. Número total de iteraciones simuladas en la cadena para muestrear la distribución a posteriori de los parámetros evolutivos. |
| **ESS (Effective Sample Size)** | Tamaño Muestral Efectivo. Métrica evaluada en Tracer para verificar la convergencia de la MCMC ($\ge 200$ recomendado para todos los parámetros). |
| **TreeAnnotator / Tracer** | Herramientas complementarias de la suite BEAST: Tracer analiza los archivos de traza `.log` y TreeAnnotator resume la muestra de árboles `.trees` en un árbol de máxima clada de credibilidad (MCC). |

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

### 1.Lanzador de CD-HIT (`lanzar_cdhit.sh`)

- **Uso:** `/ruta/carpeta/grupo/VirResp-cluster/plantillas/lanzar_cdhit.sh <secuencias.fasta>`
- **Requisito:** Archivo FASTA de nucleótidos o proteínas (no requiere alineamiento previo).

#### Inspección de entrada y autodetección
El script analiza el contenido del archivo FASTA para simplificar la configuración:
1. Cuenta el número total de secuencias y evalúa la longitud promedio.
2. Analiza automáticamente la proporción de nucleótidos vs. aminoácidos para preseleccionar la herramienta adecuada (`cd-hit-est` o `cd-hit`).

#### Configuración interactiva

1. **Tipo de Secuencia:**
   - **Nucleótidos (ADN/ARN):** Utiliza `cd-hit-est`.
   - **Aminoácidos (Proteínas):** Utiliza `cd-hit`.

2. **Umbral de Similitud de Identidad (`-c`):**
   - **100% Identidad (`-c 1.0`):** Desduplicación estricta y obtención de haplotipos únicos (ideal como paso previo a filogenia con IQ-TREE o BEAST).
   - **99% Identidad (`-c 0.99`):** Colapsado de microvariaciones y errores de PCR o secuenciación.
   - **95% Identidad (`-c 0.95`):** Agrupamiento a nivel de linajes o variantes genómicas.
   - **Personalizado:** Permite definir cualquier umbral personalizado entre `0.80` y `1.0`.

3. **Cobertura Mínima (`-aS`):**
   - Garantiza que la secuencia más corta quede alineada en el porcentaje indicado (por defecto $100\%$).

#### Estructura del Directorio de Resultados
Los resultados se organizan automáticamente en la carpeta `<basenombre>_haplotypes/`:
- `<basenombre>_representatives.fasta`: Archivo FASTA desduplicado conteniendo únicamente las secuencias representativas de cada clúster/haplotipo.
- `<basenombre>_haplotypes.txt`: Reporte legible por humanos que detalla qué muestras e identificadores específicos han sido agrupados dentro de cada haplotipo.
- `<basenombre>_raw.clstr`: Archivo nativo de agrupación de CD-HIT.
- `<basenombre>_cdhit_info.txt`: Fichero de trazabilidad con la recuento exacto de secuencias de entrada/salida, comandos aplicados y hash MD5 del archivo original.

---

### 2. Lanzador de MAFFT (`lanzar_mafft.sh`)

- **Uso:** `/ruta/carpeta/grupo/VirResp-cluster/plantillas/lanzar_mafft.sh <secuencias.fasta>`
- **Requisito:** Archivo de secuencias sin alinear en formato FASTA (*.fasta, *.fa, *.fas).

#### Inspección de entrada y validación automática
El script analiza el archivo FASTA antes de solicitar parámetros:
1. Cuenta el número total de secuencias y evalúa el rango de longitudes (mínima, máxima y promedio).
2. Detecta si el archivo ya posee longitudes idénticas y advierte al usuario si ya está alineado.
3. Detecta identificadores duplicados en las cabeceras para prevenir errores de nombrado en MAFFT.

#### Opciones de configuración interactiva

1. **Estrategia de Alineamiento:**
   - **AUTO:** Selecciona automáticamente el mejor método (`FFT-NS-1`, `FFT-NS-2`, `FFT-NS-i` o `L-INS-i`) según el tamaño del archivo.
   - **L-INS-i / G-INS-i / E-INS-i:** Estrategias iterativas de máxima precisión para conjuntos pequeños ($< 200$ secuencias).
   - **FFT-NS-2:** Método progresivo rápido para miles de secuencias.
   - **PartTree:** Algoritmo ultra-masivo diseñado para $> 10.000$ genomas completos.

2. **Ajuste de Orientación de Hebras (`--adjustdirection`):**
   - Detecta y reorienta automáticamente secuencias que se encuentren en hebra reversa complementaria (muy habitual en datos de secuenciación viral).

3. **Tratamiento de Regiones Hipervariables e Inserciones (`Leave gappy regions` / `Unalignlevel`):**
   - **Utilidad en virología:** En genomas virales con inserciones únicas (ej. inserciones en gen S de coronavirus, bucles hipervariables en genomas segmentados) o extremos $5'$ / $3'$ incompletos, los alineadores globales intentan alinear forzosamente caracteres no homólogos, distorsionando el alineamiento contiguo.
   - **Opción Activada (`--leavegappyregion --unalignlevel 0.8`):** Conserva los bloques homólogos bien conservados y deja las inserciones y zonas hipervariables sin alinear (como vacíos/gaps), evitando distorsiones filogenéticas.
   - **Nota de rendimiento:** Esta función requiere la evaluación de pares globales (`--globalpair`) y está restringida a conjuntos de hasta $1.000$ secuencias por coste de cómputo $O(N^2)$.

4. **Reordenación por Similitud (`--reorder`):**
   - Agrupa las secuencias en el archivo final por su grado de similitud genómica en lugar de mantener el orden de entrada.

5. **Penalización de Apertura de Gaps (`--op`):**
   - Permite ajustar la penalización por abrir un gap (valor por defecto estándar: 1.53).

#### Archivos Generados
Los resultados se depositan en el mismo directorio del FASTA original:
- `<basenombre>_aligned.fasta`: Alineamiento múltiple resultante en formato FASTA.
- `<basenombre>_mafft_info.txt`: Registro de reproducibilidad con la versión exacta de MAFFT, parámetros de la línea de comandos, recuento de secuencias y hash MD5 del archivo original.

---
### 3. Lanzador de FastTree (`lanzar_fasttree.sh`)

- **Uso:** `/ruta/carpeta/grupo/VirResp-cluster/plantillas/lanzar_fasttree.sh <alineamiento.fasta>`
- **Requisito:** Alineamiento en formato FASTA (*.fasta, *.aln, *.fa).

#### Comprobaciones de integridad del alineamiento
Antes de configurar el análisis, el script verifica automáticamente:
1. Longitud uniforme de las secuencias (garantiza que el alineamiento no está corrupto).
2. Ausencia de identificadores duplicados.
3. Caracteres incompatibles con formato Newick `( ) , : ; [ ]` en los nombres de las muestras.
4. Discriminación automática de tipo de datos (nucleótidos vs. aminoácidos) basada en la proporción de residuos.

#### Opciones configurables en el menú
El script permite ajustar los parámetros habituales de reconstrucción filogenética de manera explícita:

1. **Tipo de Secuencia:** `Nucleótidos` (`-nt`) o `Aminoácidos`.
2. **Modelo de Sustitución:**
   - **Nucleótidos:** `GTR+CAT` (`-gtr`) [Recomendado] o `JC+CAT` (Jukes-Cantor).
   - **Proteínas:** `LG` (`-lg`), `WAG` (`-wag`) o `JTT` (Jones-Taylor-Thornton).
3. **Categorías de Tasa de Evolución (`-cat N`):**
   - Permite seleccionar 20 (estándar), 12, 8 o fijar un valor manual (4–50).
4. **Optimización bajo Distribución Gamma (`-gamma`):**
   - `Activado (-gamma)`: Reoptimiza las longitudes de rama y calcula la verosimilitud logarítmica ($\ln L$) final bajo el modelo Gamma.
   - `Desactivado`: Mantiene las tasas estimadas directamente por la aproximación CAT.
5. **Profundidad de Búsqueda Topológica:**
   - **Estándar:** Intercambios NNI y SPR rápidos optimizados para alineamientos masivos.
   - **Exhaustiva (`-spr 4 -mlacc 2 -slownni`):** Aumenta el radio de búsqueda SPR, realiza más iteraciones NNI y aplica mayor precisión en la optimización por máxima verosimilitud.
6. **Soporte de Ramas:**
   - **SH-like Local Supports:** Evaluación rápida Shimodaira-Hasegawa en escala 0–1 (sin coste de tiempo adicional).
   - **Bootstrap Clásico (1000 repeticiones, `-boot 1000`):** Realiza 1000 pseudorréplicas de bootstrap (disponible para conjuntos $\le 5.000$ secuencias).
   - **Sin soporte (`-nosupport`).**

#### Reproducibilidad y Archivos Generados
Cada análisis fija la semilla aleatoria (`-seed 1253`) para garantizar la exactitud entre ejecuciones idénticas.

Los resultados se almacenan en la carpeta `<basenombre>_fasttree_results/`:
- `*.treefile`: Árbol filogenético en formato Newick.
- `*_parametros.txt`: Registro exhaustivo con los flags exactos empleados, versión del ejecutable, hash MD5 del archivo de entrada y cita bibliográfica.
- `*.log`: Salida estándar detallada del proceso de optimización ML.

---

### 4. Lanzador de IQ-TREE 2 (`lanzar_iqtree.sh`)

- **Uso:** `/ruta/carpeta/grupo/VirResp-cluster/plantillas/lanzar_iqtree.sh <alineamiento.fasta>`
- **Requisito:** Alineamiento múltiple en formato FASTA (*.fasta, *.aln, *.fa, *.fas).

#### Inspección del alineamiento e integridad
Antes de configurar la reconstrucción, el script ejecuta un análisis preliminar del archivo:
1. Comprueba que el archivo cuenta con un mínimo de 4 secuencias (necesarias para generar una topología filogenética no trivial).
2. Verifica estrictamente que todas las secuencias tengan la misma longitud exactas (garantiza que el archivo se encuentra alineado).
3. Inspecciona la presencia de identificadores duplicados y advierte al usuario de nombres con caracteres incompatibles con el formato Newick `( ) , : ; [ ]`.
4. Determina automáticamente la composición de residuos (nucleótidos vs. aminoácidos).

#### Configuración interactiva de parámetros

1. **Modelo de Sustitución Evolutiva:**
   - **ModelFinder Plus (`-m MFP`):** Selección automatizada del modelo evolutivo óptimo evaluando modelos Gamma y FreeRate (`+R`). [Recomendado]
   - **ModelFinder Estándar (`-m MF`):** Igual que MFP pero restringido a modelos clásicos discretos Gamma.
   - **Modelos Explícitos ADN/ARN:** `GTR+F+I+G4`, `GTR+G`, `HKY85+G`, `TN93+G`.
   - **Modelos Explícitos Aminoácidos:** `LG+F+G4`, `WAG+F+G4`, `JTT+F+G4`.
   - **Sintaxis Personalizada:** Introducción manual de cualquier modelo soportado por IQ-TREE 2.

2. **Evaluación de Soporte de Ramas (Bootstrap):**
   - **Ultrafast Bootstrap con BNNI (`-bb 1000 -bnni`):** Remuestreo ultra-rápido de 1000 réplicas con optimización NNI para evitar sobreestimación de soporte.
   - **Doble Validación (`-alrt 1000 -bb 1000 -bnni`):** Combina el test SH-aLRT con Ultrafast Bootstrap para máxima rigurosidad filogenética.
   - **Test SH-aLRT (`-alrt 1000`):** Cálculo rápido de soporte de ramas locales mediante verosimilitud aproximada.
   - **Bootstrap Estándar (`-b 100`):** Bootstrap clásico no paramétrico (reservado para alineamientos pequeños $< 500$ muestras).
   - **Sin soporte:** Estimación directa de la topología de Máxima Verosimilitud.

3. **Reproducibilidad:**
   - Permite fijar opcionalmente una semilla pseudoaleatoria fija (`-seed 1253`) para garantizar la invarianza topológica entre ejecuciones distintas.

#### Estructura del Directorio de Resultados
Cada análisis deposita sus resultados en la carpeta `<basenombre>_iqtree_results/`:
- `<basenombre>.treefile`: Árbol de Máxima Verosimilitud final en formato Newick.
- `<basenombre>.contree`: Árbol de consenso con valores de soporte de ramas en los nodos internos (cuando se activa Bootstrap/SH-aLRT).
- `<basenombre>.iqtree`: Reporte extenso con la evaluación del modelo, log-likelihood ($\ln L$), frecuencias de nucleótidos/aminoácidos y estadísticas del árbol.
- `<basenombre>.log`: Registro completo del proceso de optimización numérica.
- `<basenombre>.mldist`: Matriz de distancias evolutivas estimadas por Máxima Verosimilitud.
- `<basenombre>_iqtree_info.txt`: Fichero de auditoría con la versión ejecutable, comando exacto utilizado, hash MD5 de la entrada y marca de tiempo.

---

### 5. Lanzador de BEAST + BEAGLE (`lanzar_beast.sh`)

- **Uso:** `/ruta/carpeta/grupo/VirResp-cluster/plantillas/lanzar_beast.sh <analisis.xml>`
- **Requisito:** Archivo de configuración MCMC en formato XML generado previamente en BEAUti.

#### Inspección del archivo XML y automatización
Dado que la especificación completa del modelo de sustitución, el modelo de reloj molecular, los priores filodinámicos y las frecuencias de muestreo se definen directamente en la interfaz gráfica de BEAUti al exportar el archivo `.xml`, el script no solicita parámetros evolutivos interactivos innecesarios.

El lanzador automatiza la verificación previa:
1. Inspecciona la validez del marcado XML y verifica la presencia de etiquetas principales de BEAST.
2. Extrae automáticamente la longitud configurada de la cadena MCMC (`chainLength`) y el recuento total de taxones.
3. Evalúa la escala del problema para seleccionar automáticamente la cola de Slurm, el tiempo de ejecución y la asignación de memoria necesarias ($12\text{ h}$, $48\text{ h}$ o $5\text{ días}$).

#### Motor de aceleración computacional
El script aplica de manera predeterminada las banderas optimizadas para el entorno HPC:
- `-beagle -beagle_CPU -beagle_SSE -beagle_double`: Garantiza el uso de la librería de aceleración vectorial BEAGLE con instrucciones SSE y precisión doble, maximizando la estabilidad numérica en cálculos filodinámicos complejos.

#### Estructura del Directorio de Resultados
Los resultados de la simulación MCMC se organizan dentro de la carpeta `<basenombre>_beast_results/`:
- `<basenombre>.log`: Archivo de traza continuo tabulado con el historial de parámetros a posteriori, priors y verosimilitud (compatible con Tracer).
- `<basenombre>.trees`: Muestra de árboles filogenéticos guardados a lo largo de la MCMC (compatible con TreeAnnotator).
- `<basenombre>.ops`: Reporte con el rendimiento y grado de aceptación de las operaciones de la cadena MCMC.
- `<basenombre>_beast_info.txt`: Registro de reproducibilidad con la versión exacta de BEAST, comando de ejecución, hash MD5 del archivo `.xml` y sello temporal.n `long_idx`, tiempo máx: 5 días, 32 CPUs, 128 GB RAM.

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
