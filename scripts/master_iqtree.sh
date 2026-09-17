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

# Parseo de BASENAME limpio
BASENAME=$(basename "$INPUT_FASTA")
BASENAME="${BASENAME%.gz}"
BASENAME="${BASENAME%.fasta}"
BASENAME="${BASENAME%.fas}"
BASENAME="${BASENAME%.fa}"
BASENAME="${BASENAME%.aligned}"

# 2. CONFIGURACIÓN DEL DISCO SSD LOCAL (/local_scratch)
WORK_DIR="/local_scratch/${SLURM_JOB_ID}"
mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

# Copia usando directamente la ruta absoluta
cp "$INPUT_FASTA" input.fasta

echo "=========================================================="
echo " EJECUCIÓN EN NODO DE CÓMPUTO"
echo "=========================================================="
echo " -> Nodo asignado: $(hostname)"
echo " -> Hilos asignados: ${SLURM_CPUS_PER_TASK}"
echo " -> Comando exacto:"
echo "    iqtree2 -s input.fasta -nt AUTO -ntmax ${SLURM_CPUS_PER_TASK} ${IQTREE_EXTRA} --prefix result_tree"
echo "=========================================================="

# 3. EJECUCIÓN Y CONTROL RIGUROSO DE ERRORES
if ! iqtree2 -s input.fasta -nt AUTO -ntmax "${SLURM_CPUS_PER_TASK}" ${IQTREE_EXTRA} --prefix result_tree; then
    echo " ERROR CRÍTICO: IQ-TREE 2 ha fallado durante la reconstrucción filogenética." >&2
    echo "   Revisa el log de Slurm iqtree_${SLURM_JOB_ID}.err para conocer la causa." >&2
    
    # Intentar salvar el archivo .log generado por IQ-TREE en scratch si existe
    DEST_DIR_ERR="$(dirname "$INPUT_FASTA")/${BASENAME}_iqtree_results"
    mkdir -p "$DEST_DIR_ERR"
    [ -f result_tree.log ] && cp result_tree.log "$DEST_DIR_ERR/${BASENAME}_failed.log" || true
    
    cd /
    rm -rf "$WORK_DIR"
    exit 1
fi

# 4. RECOLECCIÓN Y RENOMBRADO ROBUSTO DE SALIDAS CRÍTICAS
DEST_DIR="$(dirname "$INPUT_FASTA")/${BASENAME}_iqtree_results"
mkdir -p "$DEST_DIR"

echo " Copiando archivos de salida generados a: $DEST_DIR"

# Mapeo explicito de archivos de salida de IQ-TREE
[ -f result_tree.treefile ] && cp result_tree.treefile "$DEST_DIR/${BASENAME}.treefile"
[ -f result_tree.contree ]  && cp result_tree.contree  "$DEST_DIR/${BASENAME}.contree"
[ -f result_tree.iqtree ]   && cp result_tree.iqtree   "$DEST_DIR/${BASENAME}.iqtree"
[ -f result_tree.log ]      && cp result_tree.log      "$DEST_DIR/${BASENAME}.log"
[ -f result_tree.mldist ]   && cp result_tree.mldist   "$DEST_DIR/${BASENAME}.mldist"

# 5. LIMPIEZA FINAL
cd /
rm -rf "$WORK_DIR"

echo "=========================================================="
echo " Análisis con IQ-TREE finalizado con éxito."
echo " Resultados consolidados en: $DEST_DIR"
echo "   - Árbol ML principal: ${BASENAME}.treefile"
echo "   - Árbol de Consenso (Bootstrap): ${BASENAME}.contree"
echo "   - Reporte detallado/Modelos: ${BASENAME}.iqtree"
echo "=========================================================="