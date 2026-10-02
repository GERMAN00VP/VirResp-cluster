#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1

# ==============================================================================
# SCRIPT MAESTRO DE EJECUCIÓN DE MAFFT EN NODO DE CÓMPUTO
# Uso: sbatch master_mafft.sh <INPUT_FASTA> "<MAFFT_EXTRA>"
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
MAFFT_EXTRA="${2:-}"

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

# Cargar entorno bash de conda de forma limpia
CONDA_BASE=$(conda info --base 2>/dev/null || echo "/opt/conda")
if [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
    source "$CONDA_BASE/etc/profile.d/conda.sh"
else
    echo " [ERROR CRÍTICO] No se pudo localizar conda.sh en $CONDA_BASE." >&2
    exit 1
fi

conda activate "$ENV_PATH"

if ! command -v mafft >/dev/null 2>&1; then
    echo " [ERROR CRÍTICO] El ejecutable 'mafft' no está disponible en el entorno $ENV_PATH." >&2
    exit 1
fi

# ------------------------------------------------------------------------------
# 3. PREPARACIÓN DE RUTAS Y ESPACIO SCRATCH LOCAL
# ------------------------------------------------------------------------------
BASENAME="$(basename "$INPUT_FASTA")"
BASENAME="${BASENAME%.fasta}"; BASENAME="${BASENAME%.fa}"; BASENAME="${BASENAME%.fas}"

ORIG_DIR="$(dirname "$(realpath "$INPUT_FASTA")")"
OUTPUT_FINAL="${ORIG_DIR}/${BASENAME}_aligned.fasta"
INFO_FILE="${ORIG_DIR}/${BASENAME}_mafft_info.txt"

# Crear directorio temporal en local scratch del nodo
WORK_DIR="/local_scratch/${SLURM_JOB_ID}"
mkdir -p "$WORK_DIR"

# Trampa de limpieza automática del directorio scratch al finalizar o abortar
trap 'rm -rf "$WORK_DIR"' EXIT

cd "$WORK_DIR"
cp "$INPUT_FASTA" input.fasta

# ------------------------------------------------------------------------------
# 4. EJECUCIÓN DEL ALINEAMIENTO
# ------------------------------------------------------------------------------
echo "=========================================================="
echo " EJECUCIÓN DE MAFFT EN NODO DE CÓMPUTO"
echo "=========================================================="
echo " -> Job ID Slurm:       ${SLURM_JOB_ID}"
echo " -> Nodo asignado:      $(hostname)"
echo " -> Hilos asignados:    ${SLURM_CPUS_PER_TASK}"
echo " -> Archivo entrada:    $INPUT_FASTA"
echo " -> Archivo salida:     $OUTPUT_FINAL"
echo " -> Banderas aplicadas: $MAFFT_EXTRA"
echo "=========================================================="

# Comando de ejecución
# Se desactiva set -e temporalmente para capturar códigos de error de MAFFT explícitamente
set +e
mafft --thread "${SLURM_CPUS_PER_TASK}" $MAFFT_EXTRA input.fasta > output.fasta 2> mafft_run.log
MAFFT_EXIT_CODE=$?
set -e

if [ $MAFFT_EXIT_CODE -ne 0 ]; then
    echo " [ERROR CRÍTICO] MAFFT falló con código de salida $MAFFT_EXIT_CODE." >&2
    echo " --- Últimas líneas del log de MAFFT ---" >&2
    tail -n 25 mafft_run.log >&2
    exit $MAFFT_EXIT_CODE
fi

# ------------------------------------------------------------------------------
# 5. VALIDACIÓN Y CONSOLIDACIÓN DE RESULTADOS
# ------------------------------------------------------------------------------
if [ ! -s output.fasta ]; then
    echo " [ERROR CRÍTICO] El alineamiento generado por MAFFT está vacío." >&2
    exit 1
fi

NUM_IN=$(grep -c '^>' input.fasta || true)
NUM_OUT=$(grep -c '^>' output.fasta || true)

if [ "$NUM_OUT" -ne "$NUM_IN" ]; then
    echo " [AVISO] El número de secuencias en la salida ($NUM_OUT) no coincide con el de entrada ($NUM_IN)." >&2
fi

# Copiar resultado final a la carpeta del usuario
cp output.fasta "$OUTPUT_FINAL"

# Generar registro de metadatos e información de reproducibilidad
{
    echo "=========================================================="
    echo " REGISTRO DE EJECUCIÓN - MAFFT ALIGNMENT"
    echo "=========================================================="
    echo " Fecha de ejecución:     $(date '+%Y-%m-%d %H:%M:%S')"
    echo " Slurm Job ID:           ${SLURM_JOB_ID}"
    echo " Archivo original:       $INPUT_FASTA"
    echo " Hash MD5 entrada:       $(md5sum "$INPUT_FASTA" | awk '{print $1}')"
    echo " Versión de MAFFT:       $(mafft --version 2>&1 || echo 'Desconocida')"
    echo " Secuencias entrada:     $NUM_IN"
    echo " Secuencias alineadas:   $NUM_OUT"
    echo " Comando exacto:         mafft --thread ${SLURM_CPUS_PER_TASK} $MAFFT_EXTRA input.fasta"
    echo "=========================================================="
} > "$INFO_FILE"

echo "=========================================================="
echo " Proceso completado con éxito."
echo " Alineamiento final:   $OUTPUT_FINAL"
echo " Registro de parámetros: $INFO_FILE"
echo "=========================================================="