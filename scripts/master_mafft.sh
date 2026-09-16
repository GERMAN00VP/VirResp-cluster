#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1

set -euo pipefail

# 1. VERIFICACIÓN DEL ENTORNO CONDA
ENV_PATH="/data/cnm/vrg/conda_envs/miniresp"

if [ ! -d "$ENV_PATH" ]; then
    echo " ERROR CRÍTICO: No se encontró el entorno Conda en '$ENV_PATH'." >&2
    exit 1
fi

source $(conda info --base)/etc/profile.d/conda.sh
conda activate "$ENV_PATH"

# Parseo robusto del BASENAME y fijación de la ruta de salida ANTES del cd
BASENAME=$(basename "$INPUT_FASTA")
BASENAME="${BASENAME%.gz}"
BASENAME="${BASENAME%.fasta}"
BASENAME="${BASENAME%.fas}"
BASENAME="${BASENAME%.fa}"

ORIG_DIR="$(dirname "$INPUT_FASTA")"
OUTPUT_FINAL="${ORIG_DIR}/${BASENAME}_aligned.fasta"

# 2. CONFIGURACIÓN DEL ESPACIO LOCAL Y TRAP DE LIMPIEZA AUTOMÁTICA
WORK_DIR="/local_scratch/${SLURM_JOB_ID}"
mkdir -p "$WORK_DIR"

# Garantiza que el espacio scratch se borre siempre al finalizar o si ocurre una falla
trap 'rm -rf "$WORK_DIR"' EXIT

cd "$WORK_DIR"
cp "$INPUT_FASTA" input.fasta

echo "=========================================================="
echo " EJECUCIÓN EN NODO DE CÓMPUTO"
echo "=========================================================="
echo " -> Nodo asignado: $(hostname)"
echo " -> Hilos asignados: ${SLURM_CPUS_PER_TASK}"
echo " -> Archivo de salida: $OUTPUT_FINAL"
echo " -> Comando exacto:"
echo "    mafft --thread ${SLURM_CPUS_PER_TASK} ${MAFFT_EXTRA} input.fasta > output.fasta"
echo "=========================================================="

# 3. EJECUCIÓN Y VALIDACIÓN ESTRICTA DEL RESULTADO
if ! mafft --thread "${SLURM_CPUS_PER_TASK}" ${MAFFT_EXTRA} input.fasta > output.fasta; then
    echo " ERROR CRÍTICO: El proceso de MAFFT ha fallado." >&2
    echo " Revisa el archivo de error mafft_${SLURM_JOB_ID}.err" >&2
    exit 1
fi

# Comprobar que el archivo resultante no esté vacío
if [ ! -s output.fasta ]; then
    echo " ERROR CRÍTICO: MAFFT ha finalizado pero ha generado un alineamiento vacío." >&2
    exit 1
fi

# Validación del número de secuencias finales
NUM_OUT=$(grep -c '^>' output.fasta || true)
if [ -n "${NUM_SEQS:-}" ] && [ "$NUM_OUT" -ne "$NUM_SEQS" ]; then
    echo " AVISO: El alineamiento generado contiene $NUM_OUT secuencias de las $NUM_SEQS del input original." >&2
fi

# 4. CONSOLIDACIÓN DE RESULTADOS
cp output.fasta "$OUTPUT_FINAL"

echo "=========================================================="
echo " Proceso finalizado correctamente."
echo " Alineamiento generado en: $OUTPUT_FINAL"
echo "=========================================================="