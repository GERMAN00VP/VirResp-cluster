#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1

# ==============================================================================
# SCRIPT MAESTRO DE EJECUCIÓN DE BEAST EN NODO DE CÓMPUTO
# Uso: sbatch master_beast.sh <INPUT_XML> "<BEAST_EXTRA>"
# ==============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# 1. COMPROBACIÓN DE ARGUMENTOS DE ENTRADA
# ------------------------------------------------------------------------------
if [ "$#" -lt 1 ]; then
    echo " [ERROR] Falta el archivo XML de entrada." >&2
    exit 1
fi

INPUT_XML="$1"
BEAST_EXTRA="${2:-}"

if [ ! -f "$INPUT_XML" ]; then
    echo " [ERROR] El archivo de entrada '$INPUT_XML' no existe." >&2
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

if ! command -v beast >/dev/null 2>&1; then
    echo " [ERROR CRÍTICO] El ejecutable 'beast' no está disponible en el entorno $ENV_PATH." >&2
    exit 1
fi

# ------------------------------------------------------------------------------
# 3. PREPARACIÓN DE RUTAS Y ESPACIO SCRATCH LOCAL
# ------------------------------------------------------------------------------
BASENAME="$(basename "$INPUT_XML" .xml)"

ORIG_DIR="$(dirname "$(realpath "$INPUT_XML")")"
DEST_DIR="${ORIG_DIR}/${BASENAME}_beast_results"
INFO_FILE="${DEST_DIR}/${BASENAME}_beast_info.txt"

WORK_DIR="/local_scratch/${SLURM_JOB_ID}"
mkdir -p "$WORK_DIR"
mkdir -p "$DEST_DIR"

# Trampa de limpieza automática del directorio scratch al finalizar o abortar
trap 'rm -rf "$WORK_DIR"' EXIT

cd "$WORK_DIR"
cp "$INPUT_XML" input.xml

# ------------------------------------------------------------------------------
# 4. EJECUCIÓN DE LA MCMC CON BEAST
# ------------------------------------------------------------------------------
echo "=========================================================="
echo " EJECUCIÓN DE BEAST EN NODO DE CÓMPUTO"
echo "=========================================================="
echo " -> Job ID Slurm:       ${SLURM_JOB_ID}"
echo " -> Nodo asignado:      $(hostname)"
echo " -> Hilos asignados:    ${SLURM_CPUS_PER_TASK}"
echo " -> Archivo entrada:    $INPUT_XML"
echo " -> Carpeta resultados: $DEST_DIR"
echo " -> Banderas aplicadas: $BEAST_EXTRA"
echo "=========================================================="

set +e
beast -threads "${SLURM_CPUS_PER_TASK}" $BEAST_EXTRA input.xml > beast_run.log 2>&1
BEAST_EXIT_CODE=$?
set -e

if [ $BEAST_EXIT_CODE -ne 0 ]; then
    echo " [ERROR CRÍTICO] BEAST falló con código de salida $BEAST_EXIT_CODE." >&2
    if [ -f beast_run.log ]; then
        echo " --- Últimas líneas del log de BEAST ---" >&2
        tail -n 30 beast_run.log >&2
        cp beast_run.log "$DEST_DIR/${BASENAME}_error.log" || true
    fi
    exit $BEAST_EXIT_CODE
fi

# ------------------------------------------------------------------------------
# 5. CONSOLIDACIÓN Y ORGANIZACIÓN DE RESULTADOS
# ------------------------------------------------------------------------------
# Recolección limpia de trazas (.log) y árboles (.trees)
if [ -f input.log ]; then
    cp input.log "$DEST_DIR/${BASENAME}.log"
else
    cp *.log "$DEST_DIR/" 2>/dev/null || true
fi

if [ -f input.trees ]; then
    cp input.trees "$DEST_DIR/${BASENAME}.trees"
else
    cp *.trees "$DEST_DIR/" 2>/dev/null || true
fi

cp *.ops "$DEST_DIR/" 2>/dev/null || true
cp *.xml.state "$DEST_DIR/" 2>/dev/null || true

# Registro de metadatos de ejecución
{
    echo "=========================================================="
    echo " REGISTRO DE EJECUCIÓN - BEAST BAYESIAN MCMC"
    echo "=========================================================="
    echo " Fecha de ejecución:     $(date '+%Y-%m-%d %H:%M:%S')"
    echo " Slurm Job ID:           ${SLURM_JOB_ID}"
    echo " Archivo original:       $INPUT_XML"
    echo " Hash MD5 entrada:       $(md5sum "$INPUT_XML" | awk '{print $1}')"
    echo " Versión ejecutable:     $(beast -version 2>&1 | head -n 1 || echo 'Desconocida')"
    echo " Comando exacto:         beast -threads ${SLURM_CPUS_PER_TASK} $BEAST_EXTRA input.xml"
    echo "=========================================================="
} > "$INFO_FILE"

echo "=========================================================="
echo " Análisis bayesiano MCMC completado con éxito."
echo " Resultados consolidados en: $DEST_DIR"
echo "   - Traza para Tracer:       ${BASENAME}.log"
echo "   - Muestra de árboles:      ${BASENAME}.trees"
echo "   - Registro de parámetros:  $INFO_FILE"
echo "=========================================================="