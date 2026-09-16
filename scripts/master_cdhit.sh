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

# Parseo de nombres de archivo y carpeta final
BASENAME=$(basename "$INPUT_FASTA")
BASENAME="${BASENAME%.gz}"
BASENAME="${BASENAME%.fasta}"
BASENAME="${BASENAME%.fas}"
BASENAME="${BASENAME%.fa}"
BASENAME="${BASENAME%.aligned}"

ORIG_DIR="$(dirname "$INPUT_FASTA")"
DEST_DIR="${ORIG_DIR}/${BASENAME}_haplotypes"
mkdir -p "$DEST_DIR"

# 2. ESPACIO SCRATCH LOCAL Y TRAP DE LIMPIEZA
WORK_DIR="/local_scratch/${SLURM_JOB_ID}"
mkdir -p "$WORK_DIR"
trap 'rm -rf "$WORK_DIR"' EXIT

cd "$WORK_DIR"
cp "$INPUT_FASTA" input.fasta

# 3. CÁLCULO ROBUSTO DE MEMORIA Y HILOS
# Cálculo compatible con SLURM_MEM_PER_NODE y SLURM_MEM_PER_CPU
RAW_MEM="${SLURM_MEM_PER_NODE:-0}"
if [ "$RAW_MEM" -eq 0 ]; then
    RAW_MEM=$(( ${SLURM_MEM_PER_CPU:-0} * ${SLURM_CPUS_PER_TASK:-1} ))
fi

if [ "$RAW_MEM" -le 0 ]; then
    RAW_MEM=32000
fi

# Margen de seguridad del 80% para evitar caídas por OOM
MEM_MB=$(( RAW_MEM * 80 / 100 ))

# Limitar hilos a máximo 8 (CD-HIT satura rendimiento con más de 8 cores)
THREADS=$(( SLURM_CPUS_PER_TASK > 8 ? 8 : SLURM_CPUS_PER_TASK ))

echo "=========================================================="
echo " EJECUCIÓN EN NODO DE CÓMPUTO"
echo "=========================================================="
echo " -> Nodo asignado: $(hostname)"
echo " -> Hilos asignados a CD-HIT: $THREADS (de ${SLURM_CPUS_PER_TASK} allocados)"
echo " -> Límite de Memoria RAM: ${MEM_MB} MB"
echo " -> Comando exacto:"
echo "    cd-hit -i input.fasta -o output_cdhit -T ${THREADS} -M ${MEM_MB} ${CDHIT_EXTRA}"
echo "=========================================================="

# 4. EJECUCIÓN Y VALIDACIÓN DE ERRORES
if ! cd-hit -i input.fasta -o output_cdhit -T "$THREADS" -M "$MEM_MB" ${CDHIT_EXTRA}; then
    echo " ERROR CRÍTICO: El proceso de CD-HIT ha fallado." >&2
    echo " Revisa el log de error cdhit_${SLURM_JOB_ID}.err" >&2
    exit 1
fi

if [ ! -s output_cdhit ]; then
    echo " ERROR CRÍTICO: CD-HIT no ha generado secuencias de salida." >&2
    exit 1
fi

# 5. GENERACIÓN DEL INFORME DETALLADO DE AGRUPACIÓN (MAPEO DE HAPLOTIPOS)
echo " Generando informe legible de haplotipos..."

awk '
/^>Cluster/ {
    if (rep != "") {
        print "HAPLOTIPO REPRESENTATIVO: " rep
        print "Total secuencias colapsadas en este grupo: " count
        print "Miembros idénticos / agrupados:"
        print seqs
        print "--------------------------------------------------------------------------------"
    }
    rep=""
    seqs=""
    count=0
    next
}
{
    # Extraer el ID exacto omitiendo > y respetando puntos internos
    match($0, />[^ ]+/)
    seq_id = substr($0, RSTART+1, RLENGTH-1)
    sub(/\.\.\.$/, "", seq_id)
    
    # Extraer el porcentaje de identidad o concordancia si existe
    match($0, /at [0-9\.]+%/)
    ident = (RLENGTH > 0) ? substr($0, RSTART, RLENGTH) : "100.00% (ref)"

    if ($0 ~ /\*$/) {
        rep = seq_id
    } else {
        seqs = (seqs == "" ? "  - " seq_id " (" ident ")" : seqs "\n  - " seq_id " (" ident ")")
    }
    count++
}
END {
    if (rep != "") {
        print "HAPLOTIPO REPRESENTATIVO: " rep
        print "Total secuencias colapsadas en este grupo: " count
        print "Miembros idénticos / agrupados:"
        print (seqs != "" ? seqs : "  (Ninguno, secuencia única)")
        print "--------------------------------------------------------------------------------"
    }
}' output_cdhit.clstr > haplotypes_summary.txt

# 6. CONSOLIDACIÓN DE ARCHIVOS DE SALIDA
cp output_cdhit "$DEST_DIR/${BASENAME}_representatives.fasta"
cp output_cdhit.clstr "$DEST_DIR/${BASENAME}_raw.clstr"
cp haplotypes_summary.txt "$DEST_DIR/${BASENAME}_haplotypes.txt"

NUM_REPS=$(grep -c '^>' output_cdhit || true)

echo "=========================================================="
echo " Proceso finalizado con éxito."
echo " Resultados guardados en: $DEST_DIR"
echo "  - FASTA desduplicado:      ${BASENAME}_representatives.fasta ($NUM_REPS haplotipos)"
echo "  - Reporte de agrupamiento: ${BASENAME}_haplotypes.txt"
echo "  - Matriz cruda CD-HIT:     ${BASENAME}_raw.clstr"
echo "=========================================================="