#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1

# ==============================================================================
# SCRIPT MAESTRO DE EJECUCIÓN DE IQ-TREE 2 EN NODO DE CÓMPUTO
# Uso: sbatch master_iqtree.sh <INPUT_FASTA> "<IQTREE_EXTRA>"
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
DEST_DIR="${ORIG_DIR}/${BASENAME}_iqtree_results"
INFO_FILE="${DEST_DIR}/${BASENAME}_iqtree_info.txt"

WORK_DIR="/local_scratch/${SLURM_JOB_ID}"
mkdir -p "$WORK_DIR"
mkdir -p "$DEST_DIR"

# Trampa de limpieza automática del directorio scratch al finalizar o abortar
trap 'rm -rf "$WORK_DIR"' EXIT

cd "$WORK_DIR"
cp "$INPUT_FASTA" input.fasta

# ------------------------------------------------------------------------------
# 4. EJECUCIÓN DEL ANÁLISIS FILOGENÉTICO
# ------------------------------------------------------------------------------
echo "=========================================================="
echo " EJECUCIÓN DE IQ-TREE 2 EN NODO DE CÓMPUTO"
echo "=========================================================="
echo " -> Job ID Slurm:       ${SLURM_JOB_ID}"
echo " -> Nodo asignado:      $(hostname)"
echo " -> Hilos asignados:    ${SLURM_CPUS_PER_TASK}"
echo " -> Archivo entrada:    $INPUT_FASTA"
echo " -> Carpeta resultados: $DEST_DIR"
echo " -> Banderas aplicadas: $IQTREE_EXTRA"
echo "=========================================================="

set +e
iqtree2 -s input.fasta -nt AUTO -ntmax "${SLURM_CPUS_PER_TASK}" $IQTREE_EXTRA --prefix result_tree > iqtree_run.log 2>&1
IQTREE_EXIT_CODE=$?
set -e

if [ $IQTREE_EXIT_CODE -ne 0 ]; then
    echo " [ERROR CRÍTICO] IQ-TREE 2 falló con código de salida $IQTREE_EXIT_CODE." >&2
    if [ -f result_tree.log ]; then
        echo " --- Últimas líneas del archivo de log de IQ-TREE ---" >&2
        tail -n 30 result_tree.log >&2
        cp result_tree.log "$DEST_DIR/${BASENAME}_error.log" || true
    fi
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
[ -f result_tree.mldist ]   && cp result_tree.mldist   "$DEST_DIR/${BASENAME}.mldist"
[ -f result_tree.model.gz ] && cp result_tree.model.gz "$DEST_DIR/${BASENAME}.model.gz"

NUM_IN=$(grep -c '^>' input.fasta || true)

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
    echo " Comando exacto:         iqtree2 -s input.fasta -nt AUTO -ntmax ${SLURM_CPUS_PER_TASK} $IQTREE_EXTRA --prefix result_tree"
    echo "=========================================================="
} > "$INFO_FILE"

echo "=========================================================="
echo " Análisis filogenético finalizado correctamente."
echo " Resultados consolidados en: $DEST_DIR"
echo "   - Árbol ML principal:     ${BASENAME}.treefile"
echo "   - Árbol de Consenso:      ${BASENAME}.contree (si se usó Bootstrap)"
echo "   - Reporte de modelos/fit: ${BASENAME}.iqtree"
echo "   - Registro parámetros:    $INFO_FILE"
echo "=========================================================="