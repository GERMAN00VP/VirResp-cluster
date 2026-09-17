#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=beast_%j.out
#SBATCH --error=beast_%j.err

set -euo pipefail

# Activar el entorno Conda que contiene BEAST y las librerías BEAGLE
source $(conda info --base)/etc/profile.d/conda.sh
conda activate /data/cnm/vrg/conda_envs/miniresp

BASENAME=$(basename "$INPUT_XML" .xml)

# Directorio de ejecución temporal rápido en el disco local del nodo SSD
WORK_DIR="/local_scratch/${SLURM_JOB_ID}"
mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

# Copiar el archivo de configuración XML al espacio de trabajo SSD
cp "$INPUT_XML" input.xml

echo "Nodo de cómputo asignado: $(hostname)"
echo "Hilos de CPU asignados: ${SLURM_CPUS_PER_TASK}"
echo "Comando exacto de ejecución:"
echo "beast -threads ${SLURM_CPUS_PER_TASK} ${BEAST_EXTRA} input.xml"

# Ejecutar BEAST
beast -threads "${SLURM_CPUS_PER_TASK}" ${BEAST_EXTRA} input.xml

# Crear carpeta ordenada de resultados finales en el directorio original del usuario
DEST_DIR="$(dirname "$INPUT_XML")/${BASENAME}_beast_results"
mkdir -p "$DEST_DIR"

# Mover los archivos generados por la MCMC (.log, .trees, .ops, etc.)
cp input.log "$DEST_DIR/${BASENAME}.log" 2>/dev/null || cp *.log "$DEST_DIR/" 2>/dev/null || true
cp input.trees "$DEST_DIR/${BASENAME}.trees" 2>/dev/null || cp *.trees "$DEST_DIR/" 2>/dev/null || true
cp *.ops "$DEST_DIR/" 2>/dev/null || true

# Limpieza del disco interno SSD
rm -rf "$WORK_DIR"

echo "=========================================================="
echo "✔ Análisis bayesiano MCMC completado con éxito."
echo "Resultados consolidados en: $DEST_DIR"
echo "  - Archivo de traza para Tracer: ${BASENAME}.log"
echo "  - Árboles muestreados para TreeAnnotator: ${BASENAME}.trees"
echo "=========================================================="