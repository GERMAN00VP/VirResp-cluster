#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
# ==============================================================================
# SCRIPT MAESTRO DE FASTTREE  (lo lanza plantillas/lanzar_fasttree.sh)
# Uso interno:  sbatch master_fasttree.sh <alineamiento.fasta> "<flags de FastTree>"
#
# Variables opcionales (solo administrador):
#   VIRRESP_ENV   ruta del entorno conda (por defecto, el del grupo)
#   FASTTREE_BIN  forzar un ejecutable concreto (p. ej. FastTreeDbl, doble precisión)
# ==============================================================================

set -euo pipefail

if [ -z "${SLURM_JOB_ID:-}" ]; then
    echo "ERROR: este script debe lanzarse con sbatch (usa plantillas/lanzar_fasttree.sh)." >&2
    exit 1
fi
if [ "$#" -lt 1 ]; then
    echo "ERROR: falta el archivo de entrada. Uso: sbatch $0 <alineamiento.fasta> \"<flags>\"" >&2
    exit 1
fi

INPUT_FASTA="$1"
FASTTREE_EXTRA="${2:-}"
SECONDS=0

if [ ! -f "$INPUT_FASTA" ]; then
    echo "ERROR: no existe el archivo de entrada: $INPUT_FASTA" >&2
    exit 1
fi

# ------------------------------------------------------------------------------
# 1. ENTORNO CONDA
# ------------------------------------------------------------------------------
ENV_PATH="${VIRRESP_ENV:-/data/cnm/vrg/conda_envs/miniresp}"

if [ ! -d "$ENV_PATH" ]; then
    echo "ERROR CRÍTICO: no se encontró el entorno conda en '$ENV_PATH'." >&2
    exit 1
fi
if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR CRÍTICO: el comando 'conda' no está disponible en el nodo de cálculo." >&2
    exit 1
fi

# Los scripts de activación de conda no son compatibles con 'set -u'
set +u
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_PATH"
set -u

# Preferir FastTreeMP (usa varios núcleos); si no existe, las versiones de 1 núcleo
FASTTREE_BIN="${FASTTREE_BIN:-}"
if [ -z "$FASTTREE_BIN" ]; then
    for cand in FastTreeMP FastTree fasttree; do
        if command -v "$cand" >/dev/null 2>&1; then
            FASTTREE_BIN="$cand"
            break
        fi
    done
fi
if [ -z "$FASTTREE_BIN" ] || ! command -v "$FASTTREE_BIN" >/dev/null 2>&1; then
    echo "ERROR CRÍTICO: no se encontró FastTree en el entorno '$ENV_PATH'." >&2
    exit 1
fi

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
if [[ "$FASTTREE_BIN" == *MP* ]]; then
    THREADS_USED="$OMP_NUM_THREADS"
else
    THREADS_USED="1 (esta versión de FastTree no usa varios núcleos)"
fi

FT_VERSION="$(conda list -p "$ENV_PATH" '^fasttree$' 2>/dev/null | awk '!/^#/ {print $1" "$2}' || true)"
FT_VERSION="${FT_VERSION:-desconocida}"

# ------------------------------------------------------------------------------
# 2. NOMBRES Y CARPETAS
# ------------------------------------------------------------------------------
BASENAME="$(basename "$INPUT_FASTA")"
BASENAME="${BASENAME%.fasta}"; BASENAME="${BASENAME%.fas}"
BASENAME="${BASENAME%.fa}";    BASENAME="${BASENAME%.fna}"

STAMP="$(date +%Y%m%d_%H%M%S)"
DEST_DIR="$(dirname "$INPUT_FASTA")/${BASENAME}_fasttree_results"
OUT_PREFIX="${BASENAME}_${STAMP}"

# Disco local rápido del nodo; si no existe, se usa un directorio temporal
SCRATCH_ROOT="/local_scratch"
if [ ! -d "$SCRATCH_ROOT" ] || [ ! -w "$SCRATCH_ROOT" ]; then
    SCRATCH_ROOT="${TMPDIR:-/tmp}"
fi
WORK_DIR="${SCRATCH_ROOT}/${SLURM_JOB_ID}_fasttree"
mkdir -p "$WORK_DIR"

# Limpieza garantizada: al terminar bien, con error, con scancel o por límite de tiempo
cleanup() { cd /; rm -rf "$WORK_DIR"; }
trap cleanup EXIT
trap 'exit 143' TERM INT

cd "$WORK_DIR"

# Copia de trabajo sin saltos de línea de Windows (el archivo original no se toca)
tr -d '\r' < "$INPUT_FASTA" > input.fasta
NUM_SEQS="$(grep -c '^>' input.fasta || true)"
INPUT_MD5="$(md5sum "$INPUT_FASTA" | cut -d' ' -f1)"

# Flags en un array para evitar problemas de espacios
read -r -a EXTRA_ARR <<< "$FASTTREE_EXTRA"

echo "=========================================================="
echo " EJECUCIÓN EN NODO DE CÓMPUTO"
echo "=========================================================="
echo " -> Nodo:        $(hostname)"
echo " -> Secuencias:  $NUM_SEQS"
echo " -> Programa:    $FASTTREE_BIN (paquete: $FT_VERSION)"
echo " -> Hilos:       $THREADS_USED"
echo " -> Comando:     $FASTTREE_BIN -log result.log ${FASTTREE_EXTRA} input.fasta > result.treefile"
echo "=========================================================="

# ------------------------------------------------------------------------------
# 3. EJECUCIÓN
# ------------------------------------------------------------------------------
if ! "$FASTTREE_BIN" -log result.log ${EXTRA_ARR[@]+"${EXTRA_ARR[@]}"} input.fasta > result.treefile \
   || [ ! -s result.treefile ] || ! grep -q ';' result.treefile; then
    echo "" >&2
    echo "ERROR CRÍTICO: FastTree no ha podido generar el árbol." >&2
    echo "  Mira las últimas líneas de fasttree_${SLURM_JOB_ID}.err para ver la causa." >&2
    mkdir -p "$DEST_DIR"
    if [ -f result.log ]; then
        cp result.log "$DEST_DIR/${OUT_PREFIX}_FALLO.log" || true
        echo "  Se ha guardado el registro parcial en: $DEST_DIR/${OUT_PREFIX}_FALLO.log" >&2
    fi
    exit 1
fi

# ------------------------------------------------------------------------------
# 4. REGISTRO DE PARÁMETROS (para poder reproducir y citar el análisis)
# ------------------------------------------------------------------------------
{
    echo "VirResp-cluster | Registro de análisis con FastTree"
    echo "---------------------------------------------------"
    echo "Fecha de finalización : $(date '+%Y-%m-%d %H:%M:%S')"
    echo "Duración              : ${SECONDS} s"
    echo "Trabajo Slurm (ID)    : $SLURM_JOB_ID"
    echo "Nodo                  : $(hostname)"
    echo "Archivo de entrada    : $INPUT_FASTA"
    echo "MD5 de la entrada     : $INPUT_MD5"
    echo "Nº de secuencias      : $NUM_SEQS"
    echo "Programa              : $FASTTREE_BIN"
    echo "Versión (paquete)     : $FT_VERSION"
    echo "Hilos                 : $THREADS_USED"
    echo "Parámetros            : $FASTTREE_EXTRA"
    echo "Comando               : $FASTTREE_BIN -log ${OUT_PREFIX}.log $FASTTREE_EXTRA ${BASENAME}.fasta > ${OUT_PREFIX}.treefile"
    echo "Cita                  : Price MN, Dehal PS, Arkin AP (2010) FastTree 2 - Approximately"
    echo "                        Maximum-Likelihood Trees for Large Alignments. PLoS ONE 5(3): e9490."
} > params.txt

# ------------------------------------------------------------------------------
# 5. RECOGIDA DE RESULTADOS (sin sobrescribir ejecuciones anteriores)
# ------------------------------------------------------------------------------
mkdir -p "$DEST_DIR"
cp result.treefile "$DEST_DIR/${OUT_PREFIX}.treefile"
cp params.txt      "$DEST_DIR/${OUT_PREFIX}_parametros.txt"
if [ -f result.log ]; then
    cp result.log "$DEST_DIR/${OUT_PREFIX}.log"
fi

if [ ! -s "$DEST_DIR/${OUT_PREFIX}.treefile" ]; then
    echo "ERROR CRÍTICO: no se pudo copiar el árbol a '$DEST_DIR'." >&2
    exit 1
fi

echo "=========================================================="
echo " Análisis con FastTree finalizado con éxito."
echo " Resultados en: $DEST_DIR"
echo "   - Árbol (sin raíz):     ${OUT_PREFIX}.treefile"
echo "   - Parámetros usados:    ${OUT_PREFIX}_parametros.txt"
echo "   - Registro detallado:   ${OUT_PREFIX}.log"
echo " Los números sobre las ramas van de 0 a 1 (0,95 = 95 %)."
echo "=========================================================="
