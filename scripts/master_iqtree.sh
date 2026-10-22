#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1

# ==============================================================================
# SCRIPT MAESTRO DE EJECUCIÓN DE IQ-TREE 2 EN NODO DE CÓMPUTO
# Uso: sbatch master_iqtree.sh <INPUT_FASTA> "<IQTREE_EXTRA>" [OUTDIR]
# ==============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# 1. COMPROBACIÓN DE ARGUMENTOS DE ENTRADA
# ------------------------------------------------------------------------------
if [ "$#" -lt 1 ]; then
    echo " [ERROR] Falta el archivo FASTA de entrada." >&2
    exit 1
fi

INPUT_FASTA="$1"
IQTREE_EXTRA="${2:-}"
OUTDIR_ARG="${3:-}"

if [ ! -f "$INPUT_FASTA" ]; then
    echo " [ERROR] El archivo de entrada '$INPUT_FASTA' no existe." >&2
    exit 1
fi

# ------------------------------------------------------------------------------
# 2. ACTIVACIÓN DEL ENTORNO CONDA
# ------------------------------------------------------------------------------
ENV_PATH="/data/cnm/vrg/conda_envs/miniresp"

if [ ! -d "$ENV_PATH" ]; then
    echo " [ERROR CRÍTICO] No se encontró el entorno Conda en '$ENV_PATH'." >&2
    exit 1
fi

CONDA_BASE=$(conda info --base 2>/dev/null || echo "/opt/conda")
if [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
    source "$CONDA_BASE/etc/profile.d/conda.sh"
else
    echo " [ERROR CRÍTICO] No se pudo localizar conda.sh en $CONDA_BASE." >&2
    exit 1
fi

conda activate "$ENV_PATH"

if ! command -v iqtree2 >/dev/null 2>&1; then
    echo " [ERROR CRÍTICO] El ejecutable 'iqtree2' no está disponible en el entorno $ENV_PATH." >&2
    exit 1
fi

# ------------------------------------------------------------------------------
# 3. PREPARACIÓN DE RUTAS Y ESPACIO SCRATCH LOCAL
# ------------------------------------------------------------------------------
BASENAME="$(basename "$INPUT_FASTA")"
BASENAME="${BASENAME%.fasta}";  BASENAME="${BASENAME%.aln}"
BASENAME="${BASENAME%.fas}";    BASENAME="${BASENAME%.fa}"
BASENAME="${BASENAME%.fna}";    BASENAME="${BASENAME%.aligned}"

ORIG_DIR="$(dirname "$(realpath "$INPUT_FASTA")")"
# Directorio elegido en el lanzador; si no se pasa, el comportamiento antiguo
DEST_DIR="${OUTDIR_ARG:-${ORIG_DIR}/${BASENAME}_iqtree_results}"
INFO_FILE="${DEST_DIR}/${BASENAME}_iqtree_info.txt"

RUN_LOG="${DEST_DIR}/${BASENAME}_iqtree_run.log"
MAX_COPY_BYTES=1073741824   # no se copian ficheros auxiliares de más de 1 GB

WORK_DIR="/local_scratch/${SLURM_JOB_ID}"
mkdir -p "$WORK_DIR"
mkdir -p "$DEST_DIR"

# Antes de borrar el scratch se copian los resultados parciales (incluido el
# checkpoint .ckp.gz) para poder reanudar si el job se corta o falla.
COMPLETADO=0
SYNC_PID=""
salvar() {
    local f
    [ "$COMPLETADO" -eq 1 ] && return 0
    for f in "$WORK_DIR"/result_tree.*; do
        [ -f "$f" ] || continue
        if [ "$(stat -c %s "$f")" -gt "$MAX_COPY_BYTES" ]; then
            echo " [AVISO] No se copia $(basename "$f") (> 1 GB)."
            continue
        fi
        cp -f "$f" "$DEST_DIR"/ 2>/dev/null || true
    done
}
trap 'salvar; [ -n "$SYNC_PID" ] && kill "$SYNC_PID" 2>/dev/null; rm -rf "$WORK_DIR"' EXIT
trap 'echo " [AVISO] Límite de tiempo cercano: copiando resultados parciales."; salvar' USR1
trap 'exit 143' TERM INT

cd "$WORK_DIR"
# Copia de trabajo sin saltos de línea de Windows (el archivo original no se toca)
tr -d '\r' < "$INPUT_FASTA" > input.fasta

# Hilos fijos = CPUs reservadas, salvo que el lanzador ya pase -T AUTO --threads-max N
if [[ " $IQTREE_EXTRA " == *" -T "* ]]; then
    THREAD_FLAGS=""
else
    THREAD_FLAGS="-T ${SLURM_CPUS_PER_TASK}"
fi

# ------------------------------------------------------------------------------
# 4. EJECUCIÓN DEL ANÁLISIS FILOGENÉTICO
# ------------------------------------------------------------------------------
echo "=========================================================="
echo " EJECUCIÓN DE IQ-TREE 2 EN NODO DE CÓMPUTO"
echo "=========================================================="
echo " -> Log en tiempo real: $RUN_LOG"
echo "    (síguelo con: tail -f $RUN_LOG)"
echo " -> Job ID Slurm:       ${SLURM_JOB_ID}"
echo " -> Nodo asignado:      $(hostname)"
echo " -> CPUs reservadas:    ${SLURM_CPUS_PER_TASK}"
echo " -> Hilos IQ-TREE:      ${THREAD_FLAGS:-AUTO (en las banderas)}"
echo " -> Archivo entrada:    $INPUT_FASTA"
echo " -> Carpeta resultados: $DEST_DIR"
echo " -> Banderas aplicadas: $IQTREE_EXTRA"
echo "=========================================================="

# Copia periódica del checkpoint al directorio de salida (cada 30 min)
( while sleep 1800; do cp -f "$WORK_DIR"/result_tree.ckp.gz "$DEST_DIR"/ 2>/dev/null || true; done ) &
SYNC_PID=$!

# En segundo plano + wait: así el trap USR1 se atiende mientras IQ-TREE corre
set +e
iqtree2 -s input.fasta $THREAD_FLAGS $IQTREE_EXTRA --prefix result_tree \
    > "$RUN_LOG" 2>&1 &
TOOL_PID=$!
IQTREE_EXIT_CODE=0
while kill -0 "$TOOL_PID" 2>/dev/null; do
    wait "$TOOL_PID"; IQTREE_EXIT_CODE=$?
done
set -e
kill "$SYNC_PID" 2>/dev/null || true
SYNC_PID=""

if [ $IQTREE_EXIT_CODE -ne 0 ]; then
    echo " [ERROR CRÍTICO] IQ-TREE 2 falló con código de salida $IQTREE_EXIT_CODE." >&2
    if [ -f result_tree.log ]; then
        echo " --- Últimas líneas del archivo de log de IQ-TREE ---" >&2
        tail -n 30 result_tree.log >&2
        cp result_tree.log "$DEST_DIR/${BASENAME}_error.log" || true
    fi
    echo " Los ficheros parciales (result_tree.*, incluido el checkpoint) se copian a $DEST_DIR" >&2
    exit $IQTREE_EXIT_CODE
fi

# ------------------------------------------------------------------------------
# 5. CONSOLIDACIÓN Y ORGANIZACIÓN DE RESULTADOS
# ------------------------------------------------------------------------------
if [ ! -f result_tree.treefile ]; then
    echo " [ERROR CRÍTICO] IQ-TREE 2 no generó el archivo principal de árbol (.treefile)." >&2
    exit 1
fi

# Mapeo explicito y seguro de las salidas principales de IQ-TREE
[ -f result_tree.treefile ] && cp result_tree.treefile "$DEST_DIR/${BASENAME}.treefile"
[ -f result_tree.contree ]  && cp result_tree.contree  "$DEST_DIR/${BASENAME}.contree"
[ -f result_tree.iqtree ]   && cp result_tree.iqtree   "$DEST_DIR/${BASENAME}.iqtree"
[ -f result_tree.log ]      && cp result_tree.log      "$DEST_DIR/${BASENAME}.log"
[ -f result_tree.model.gz ] && cp result_tree.model.gz "$DEST_DIR/${BASENAME}.model.gz"
# La matriz de distancias puede ocupar varios GB con miles de secuencias
if [ -f result_tree.mldist ]; then
    if [ "$(stat -c %s result_tree.mldist)" -gt "$MAX_COPY_BYTES" ]; then
        echo " [AVISO] No se copia ${BASENAME}.mldist: ocupa $(du -h result_tree.mldist | cut -f1) (> 1 GB)." | tee -a "$RUN_LOG"
    else
        cp result_tree.mldist "$DEST_DIR/${BASENAME}.mldist"
    fi
fi
COMPLETADO=1

NUM_IN=$(grep -c '^>' input.fasta || true)

# Hilos reales: con -T AUTO, IQ-TREE escribe el número elegido en su log
THREADS_REAL="${SLURM_CPUS_PER_TASK}"
if [ -z "$THREAD_FLAGS" ]; then
    THREADS_REAL=$(grep -m1 -oE 'BEST NUMBER OF THREADS: [0-9]+' result_tree.log 2>/dev/null | grep -oE '[0-9]+$' || true)
    THREADS_REAL="${THREADS_REAL:-AUTO (no detectado)}"
fi

# Registro de métricas e información de reproducibilidad
{
    echo "=========================================================="
    echo " REGISTRO DE EJECUCIÓN - IQ-TREE 2 PHYLOGENY"
    echo "=========================================================="
    echo " Fecha de ejecución:     $(date '+%Y-%m-%d %H:%M:%S')"
    echo " Slurm Job ID:           ${SLURM_JOB_ID}"
    echo " Archivo original:       $INPUT_FASTA"
    echo " Hash MD5 entrada:       $(md5sum "$INPUT_FASTA" | awk '{print $1}')"
    echo " Versión ejecutable:     $(iqtree2 --version | head -n 1)"
    echo " Secuencias evaluadas:   $NUM_IN"
    echo " Hilos usados:           $THREADS_REAL (CPUs reservadas: ${SLURM_CPUS_PER_TASK})"
    echo " Directorio de salida:   $DEST_DIR"
    echo " Log de ejecución:       $RUN_LOG"
    echo " Comando exacto:         iqtree2 -s input.fasta $THREAD_FLAGS $IQTREE_EXTRA --prefix result_tree"
    echo "=========================================================="
} > "$INFO_FILE"

echo "=========================================================="
echo " Análisis filogenético finalizado correctamente."
echo " Resultados consolidados en: $DEST_DIR"
echo "   - Árbol ML principal:     ${BASENAME}.treefile"
echo "   - Árbol de Consenso:      ${BASENAME}.contree (si se usó Bootstrap)"
echo "   - Reporte de modelos/fit: ${BASENAME}.iqtree"
echo "   - Log de ejecución:       $RUN_LOG"
echo "   - Registro parámetros:    $INFO_FILE"
echo "=========================================================="
