```markdown
<!-- INICIO DE DOCUMENTO README.md -->

# Manual del Clúster HPC - Grupo de Virus Respiratorios (VirResp)

Bienvenido al repositorio central de herramientas bioinformáticas del grupo VirResp. Este sistema está diseñado para que puedas ejecutar análisis de alineamiento, desduplicación y filogenia sin necesidad de programar ni configurar manualmente el gestor de colas del clúster.

---

## Índice de Navegación Rápida

### Guía Rápida
* [Comandos exprés de Slurm](#comandos-exprés-de-slurm)
* [Estructura del Repositorio](#estructura-del-repositorio)
* [Flujo de Trabajo Interactivo en 3 Pasos](#flujo-de-trabajo-interactivo-en-3-pasos)
* [Subir y Bajar Datos del Clúster](#subir-y-bajar-datos-del-clúster)

---

### Guía Detallada
1. [Antes de Empezar: Requisitos](#1-antes-de-empezar-requisitos)
2. [Comandos Básicos de Linux](#2-comandos-básicos-de-linux)
3. [Cómo Funciona el Clúster y Slurm](#3-cómo-funciona-el-clúster-y-slurm)
4. [Entorno Conda: ¿Qué necesito activar?](#4-entorno-conda-qué-necesito-activar)
5. [Guía de Herramientas Disponibles](#5-guía-de-herramientas-disponibles)
   * [MAFFT (Alineamiento de Secuencias)](#mafft-alineamiento-de-secuencias)
   * [CD-HIT (Haplotipado y Desduplicación)](#cd-hit-haplotipado-y-desduplicación)
   * [IQ-TREE (Filogenia por Máxima Verosimilitud)](#iq-tree-filogenia-por-máxima-verosimilitud)
   * [BEAST (Relojes Moleculares y Filodinámica)](#beast-relojes-moleculares-y-filodinámica)
6. [Resolución de Errores Típicos](#6-resolución-de-errores-típicos)

---

## Comandos exprés de Slurm

`squeue -u $USER` : Muestra tus trabajos en ejecución o en espera.  
`scancel <JOB_ID>` : Cancela un trabajo específico (ejemplo: `scancel 123456`).  
`scancel -u $USER` : Cancela TODOS tus trabajos activos.  
`sinfo` : Muestra la disponibilidad de los nodos del clúster.  

---

## Estructura del Repositorio

El proyecto `VirResp-cluster` está organizado en dos niveles para evitar modificaciones accidental del código de cómputo:

<!-- INICIO BLOQUE TEXT: ESTRUCTURA ARCHIVOS -->
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
<!-- FIN BLOQUE TEXT: ESTRUCTURA ARCHIVOS -->

### Arquitectura de Ejecución

<!-- INICIO BLOQUE TEXT: DIAGRAMA ARQUITECTURA -->
```text
Usuario en Terminal
      │
      ▼
Ejecuta: plantillas/lanzar_herramienta.sh archivo.fasta
      │
      ├─► 1. Analiza el archivo (cuenta secuencias, longitud)
      ├─► 2. Muestra un menú de preguntas bioinformáticas
      ├─► 3. Configura los parámetros y los recursos de Slurm
      └─► 4. Envía el trabajo mediante 'sbatch' al Script Maestro
                  │
                  ▼
          scripts/master_herramienta.sh (en nodo de cómputo)
                  │
                  ├─► Activa el entorno Conda ('miniresp')
                  ├─► Copia los datos a /local_scratch
                  ├─► Ejecuta el análisis bioinformático
                  └─► Copia los resultados a la carpeta del usuario
```
<!-- FIN BLOQUE TEXT: DIAGRAMA ARQUITECTURA -->

---

## Flujo de Trabajo Interactivo en 3 Pasos

No necesitas editar variables dentro del código ni modificar archivos de texto para enviar trabajos al clúster.

### Paso 1: Navega a la carpeta donde están tus datos
<!-- INICIO BLOQUE BASH: PASO 1 -->
```bash
cd /data/cnm/tu_usuario/mi_analisis/
```
<!-- FIN BLOQUE BASH: PASO 1 -->

### Paso 2: Ejecuta el script interactivo correspondiente
Pasa tu archivo de entrada como argumento al script interactivo:
<!-- INICIO BLOQUE BASH: PASO 2 -->
```bash
~/Documentos/VirResp-cluster/plantillas/lanzar_mafft.sh mis_secuencias.fasta
```
<!-- FIN BLOQUE BASH: PASO 2 -->

### Paso 3: Responde a las preguntas en pantalla
El script desplegará un menú interactivo. Selecciona las opciones deseadas respondiendo con números o letras (ej. '1', 'S', 'N') y confirma el envío al clúster.

---

## Subir y Bajar Datos del Clúster

### 1. Desde una interfaz gráfica (Recomendado para uso diario)
Puedes utilizar programas como **WinSCP** (Windows) o **FileZilla** (Mac/Linux):
* **Host:** IP o dominio del clúster.
* **Puerto:** 22 (o el puerto SSH asignado).
* **Protocolo:** SFTP.
* Navega en el panel derecho hasta tu directorio personal `/data/cnm/tu_usuario/` para arrastrar y soltar archivos.

### 2. Desde la línea de comandos (Terminal local)

* **Subir un archivo al clúster:**
<!-- INICIO BLOQUE BASH: SUBIR SCP -->
```bash
scp mis_secuencias.fasta usuario@cluster:/data/cnm/tu_usuario/mi_analisis/
```
<!-- FIN BLOQUE BASH: SUBIR SCP -->

* **Descargar los resultados a tu ordenador:**
<!-- INICIO BLOQUE BASH: BAJAR SCP -->
```bash
scp -r usuario@cluster:/data/cnm/tu_usuario/mi_analisis/resultados/ ./
```
<!-- FIN BLOQUE BASH: BAJAR SCP -->

* **Sincronizar carpetas grandes (rsync):**
<!-- INICIO BLOQUE BASH: RSYNC -->
```bash
rsync -avz --progress usuario@cluster:/data/cnm/tu_usuario/mi_analisis/ ./
```
<!-- FIN BLOQUE BASH: RSYNC -->

---

## 1. Antes de Empezar: Requisitos

Antes de lanzar cualquier proceso en el clúster, asegúrate de cumplir con lo siguiente:
1. Tener una cuenta activa de usuario en el clúster.
2. Estar conectado a la red local del centro o a la VPN institucional.
3. Tener permisos de escritura en la ruta del laboratorio: `/data/cnm/vrg/` o en tu espacio `/data/cnm/tu_usuario/`.
4. Clonar este repositorio en tu directorio personal si no está presente:
<!-- INICIO BLOQUE BASH: GIT CLONE -->
```bash
git clone https://github.com/tu_repositorio/VirResp-cluster.git ~/Documentos/VirResp-cluster
```
<!-- FIN BLOQUE BASH: GIT CLONE -->

---

## 2. Comandos Básicos de Linux

| Comando | Descripción | Ejemplo de uso |
| :--- | :--- | :--- |
| `pwd` | Muestra la ruta completa del directorio actual | `pwd` |
| `ls -lh` | Lista los archivos del directorio mostrando su tamaño | `ls -lh` |
| `cd <ruta>` | Cambia al directorio especificado | `cd /data/cnm/tu_usuario/` |
| `mkdir <nombre>` | Crea una carpeta nueva | `mkdir brote_influ_2026` |
| `cp <origen> <destino>` | Copia un archivo o directorio | `cp datos.fasta ./carpeta/` |
| `grep -c "^>" <archivo>` | Cuenta cuántas secuencias FASTA existen | `grep -c "^>" secuencias.fasta` |
| `head -n 20 <archivo>` | Muestra las primeras 20 líneas de un archivo | `head -n 20 archivo.txt` |

---

## 3. Cómo Funciona el Clúster y Slurm

El clúster procesa los datos mediante un gestor de colas llamado **Slurm**. No se deben ejecutar análisis directamente en el nodo de acceso (login node) porque se bloquea el sistema para los demás usuarios.

### Ubicación de los archivos de salida y logs
Al ejecutar los scripts interactivos, los archivos de registro (`.out` y `.err`) se guardarán automáticamente **en el mismo directorio donde se encuentra tu archivo de entrada**:
* `herramienta_<JOB_ID>.out`: Muestra el registro paso a paso y el estado del análisis.
* `herramienta_<JOB_ID>.err`: Muestra alertas técnicas o mensajes de error si el proceso falla.

---

## 4. Entorno Conda: ¿Qué necesito activar?

**No necesitas instalar nada ni activar entornos de software manualmente.**

Los scripts maestros activan de forma automática el entorno unificado de Conda (`miniresp`), configurado por Germán. Este entorno contiene las versiones exactas de MAFFT, CD-HIT, IQ-TREE, BEAST y sus dependencias.

*Nota:* El archivo `environment.yml` presente en la raíz del repositorio se utiliza únicamente para tareas de mantenimiento o reconstrucción del entorno por parte del administrador.

---

## 5. Guía de Herramientas Disponibles

### MAFFT (Alineamiento de Secuencias)

Alineamiento múltiple de secuencias genómicas o de genes individuales.

* **Ejecución:**
<!-- INICIO BLOQUE BASH: EJECUCION MAFFT -->
```bash
~/Documentos/VirResp-cluster/plantillas/lanzar_mafft.sh mi_fichero.fasta
```
<!-- FIN BLOQUE BASH: EJECUCION MAFFT -->

* **Opciones del menú:**
  * **AUTO:** Opción por defecto recomendada para la mayoría de los análisis.
  * **L-INS-i:** Máxima precisión para regiones hipervariables (< 500 secuencias).
  * **FFT-NS-2:** Algoritmo rápido para conjuntos de datos grandes.
  * **PartTree:** Algoritmo masivo cuando se trabaja con más de 10.000 genomas.
* **Aspecto clave en virología:** Permite activar la opción `--adjustdirection`. Se recomienda marcar **S** si las muestras pueden contener secuencias alineadas en sentido reverso-complementario.
* **Archivos generados:**
  * `<nombre>_aligned.fasta`: Archivo FASTA alineado final.

---

### CD-HIT (Haplotipado y Desduplicación)

Identificación de secuencias idénticas o muy similares para reducir la redundancia en conjuntos de datos virológicos.

* **Ejecución:**
<!-- INICIO BLOQUE BASH: EJECUCION CDHIT -->
```bash
~/Documentos/VirResp-cluster/plantillas/lanzar_cdhit.sh mi_fichero.fasta
```
<!-- FIN BLOQUE BASH: EJECUCION CDHIT -->

* **Opciones del menú:**
  * **100% Identidad (-c 1.0):** Haplotipado estricto. Obtiene únicamente secuencias 100% idénticas.
  * **99% Identidad (-c 0.99):** Desduplicación permitiendo pequeños errores de secuenciación o PCR.
  * **95% Identidad (-c 0.95):** Agrupamiento por linajes o genotipos cercanos.
* **Archivos generados:**
  Crea la carpeta `<nombre>_haplotypes/` con los siguientes archivos:
  * `<nombre>_representatives.fasta`: Secuencias desduplicadas (un representante por grupo).
  * `<nombre>_haplotypes.txt`: Reporte detallado en texto plano que lista qué secuencia representa a qué grupo.
  * `<nombre>_raw.clstr`: Matriz de agrupación nativa de CD-HIT.

---

### IQ-TREE (Filogenia por Máxima Verosimilitud)

Construcción de árboles filogenéticos bajo el criterio de Máxima Verosimilitud.

* **Ejecución:**
<!-- INICIO BLOQUE BASH: EJECUCION IQTREE -->
```bash
~/Documentos/VirResp-cluster/plantillas/lanzar_iqtree.sh mi_alineamiento.fasta
```
<!-- FIN BLOQUE BASH: EJECUCION IQTREE -->

* **Opciones del menú:**
  * **ModelFinder:** Evaluación automática del mejor modelo de sustitución nucleotídica.
  * **UFBoot:** Cálculo del soporte de ramas mediante Ultra-Fast Bootstrap (1000 réplicas).
* **Archivos generados:**
  * `<nombre>.treefile`: Árbol filogenético final en formato Newick (compatible con FigTree, iTOL, Microreact o Auspice).
  * `<nombre>.iqtree`: Reporte detallado del modelo de sustitución seleccionado y estadísticos del árbol.
  * `<nombre>.log`: Log con el avance numérico de la optimización del árbol.

---

### BEAST (Relojes Moleculares y Filodinámica)

Análisis filogenético bayesiano para la estimación de fechas de divergencia y tasas de evolución.

* **Flujo previo fuera del clúster:**
  1. Abre la aplicación de escritorio **BEAUti** en tu ordenador personal.
  2. Carga tu alineamiento FASTA o NEXUS y configura las fechas de muestreo, el modelo de sustitución, el reloj molecular y el modelo poblacional.
  3. Exporta el archivo de configuración `.xml`.
  4. Sube **únicamente el archivo .xml** al clúster (no subas el archivo FASTA fuente a la carpeta de BEAST).

* **Ejecución en el clúster:**
<!-- INICIO BLOQUE BASH: EJECUCION BEAST -->
```bash
~/Documentos/VirResp-cluster/plantillas/lanzar_beast.sh mi_analisis.xml
```
<!-- FIN BLOQUE BASH: EJECUCION BEAST -->

* **Archivos generados:**
  * `<nombre>.log`: Tabla de parámetros muestreados durante la cadena MCMC. Se analiza en tu ordenador con la herramienta **Tracer**.
  * `<nombre>.trees`: Conjunto de árboles muestreados durante la ejecución de la cadena MCMC.

---

## 6. Resolución de Errores Típicos

### ❌ CANCELLED DUE TO TIME LIMIT
* **Causa:** El tiempo de ejecución superó el máximo asignado a la cola de Slurm.
* **Solución:** Los scripts ajustan los tiempos automáticamente según el volumen de datos. Si tu trabajo vence por tiempo, selecciona algoritmos más rápidos en el menú interactivo (ejemplo: usar `FFT-NS-2` en lugar de `L-INS-i` en MAFFT) o contacta con Germán para solicitar una extensión a la partición de larga duración (`long_idx`).

### ❌ Out Of Memory / OOM Killed
* **Causa:** El proceso ha superado la memoria RAM reservada.
* **Solución:** Los scripts ajustan los límites de memoria dinámicamente según la herramienta (ejemplo: CD-HIT detecta conjuntos de datos de > 50.000 secuencias y MAFFT a partir de 10.000). Si ocurre este fallo, desduplica primero el conjunto de muestras usando CD-HIT al 100% de identidad antes de procesarlo en MAFFT o IQ-TREE.

### ❌ Caracteres invisibles de Windows (`\r`)
* **Causa:** El archivo de datos se creó o editó en un sistema Windows y conserva saltos de línea incompatibles con Linux.
* **Solución:** Los scripts interactivos corrigen el archivo automáticamente antes de procesarlo. Si deseas solucionarlo manualmente:
<!-- INICIO BLOQUE BASH: FIX WINDOWS -->
```bash
sed -i 's/\r$//' mi_fichero.fasta
```
<!-- FIN BLOQUE BASH: FIX WINDOWS -->

### ❌ Permission denied
* **Causa:** Los scripts interactivos de la carpeta `plantillas/` no disponen de permisos de ejecución en Linux.
* **Solución:** Asigna permisos de ejecución mediante el comando:
<!-- INICIO BLOQUE BASH: FIX CHMOD -->
```bash
chmod +x ~/Documentos/VirResp-cluster/plantillas/*.sh
```
<!-- FIN BLOQUE BASH: FIX CHMOD -->

---

*Repositorio mantenido por Germán. Ante cualquier duda, incidencia o solicitud de nuevas herramientas, contacta directamente con el administrador.*

<!-- FIN DE DOCUMENTO README.md -->
```