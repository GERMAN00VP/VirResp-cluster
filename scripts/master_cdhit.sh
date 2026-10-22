#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1

# ==============================================================================
# SCRIPT MAESTRO DE EJECUCIÓN DE CD-HIT EN NODO DE CÓMPUTO
# Uso: sbatch master_cdhit.sh <INPUT_FASTA> <CDHIT_CMD> "<CDHIT_EXTRA>" [OUTDIR]
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
CDHIT_CMD="${2:-cd-hit-est}"
CDHIT_EXTRA="${3:-}"
OUTDIR_ARG="${4:-}"

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

if ! command -v "$CDHIT_CMD" >/dev/null 2>&1; then
    echo " [ERROR CRÍTICO] El ejecutable '$CDHIT_CMD' no está disponible en el entorno $ENV_PATH." >&2
    exit 1
fi

# ------------------------------------------------------------------------------
# 3. PREPARACIÓN DE RUTAS Y ESPACIO SCRATCH LOCAL
# ------------------------------------------------------------------------------
BASENAME="$(basename "$INPUT_FASTA")"
BASENAME="${BASENAME%.fasta}"; BASENAME="${BASENAME%.fa}"
BASENAME="${BASENAME%.fas}";   BASENAME="${BASENAME%.aligned}"

ORIG_DIR="$(dirname "$(realpath "$INPUT_FASTA")")"
# Directorio elegido en el lanzador; si no se pasa, el comportamiento antiguo
DEST_DIR="${OUTDIR_ARG:-${ORIG_DIR}/${BASENAME}_haplotypes}"
INFO_FILE="${DEST_DIR}/${BASENAME}_cdhit_info.txt"
RUN_LOG="${DEST_DIR}/${BASENAME}_cdhit_run.log"

WORK_DIR="/local_scratch/${SLURM_JOB_ID}"
mkdir -p "$WORK_DIR"
mkdir -p "$DEST_DIR"

# Antes de borrar el scratch se copia lo que haya (resultado parcial) si el job no terminó
COMPLETADO=0
salvar() {
    [ "$COMPLETADO" -eq 1 ] && return 0
    [ -s "$WORK_DIR/output_cdhit" ] && cp -f "$WORK_DIR/output_cdhit" "$DEST_DIR/${BASENAME}_representatives.PARCIAL.fasta" 2>/dev/null || true
    [ -s "$WORK_DIR/output_cdhit.clstr" ] && cp -f "$WORK_DIR/output_cdhit.clstr" "$DEST_DIR/${BASENAME}_raw.PARCIAL.clstr" 2>/dev/null || true
    return 0
}
trap 'salvar; rm -rf "$WORK_DIR"' EXIT
trap 'echo " [AVISO] Límite de tiempo cercano: copiando resultados parciales."; salvar' USR1
trap 'exit 143' TERM INT

cd "$WORK_DIR"
# Copia de trabajo sin saltos de línea de Windows (el archivo original no se toca)
tr -d '\r' < "$INPUT_FASTA" > input.fasta

# ------------------------------------------------------------------------------
# 4. CÁLCULO OPTIMIZADO DE RECURSOS (RAM Y HILOS)
# ------------------------------------------------------------------------------
RAW_MEM="${SLURM_MEM_PER_NODE:-0}"
if [ "$RAW_MEM" -eq 0 ]; then
    RAW_MEM=$(( ${SLURM_MEM_PER_CPU:-0} * ${SLURM_CPUS_PER_TASK:-1} ))
fi

if [ "$RAW_MEM" -le 0 ]; then
    RAW_MEM=32000
fi

# Margen de seguridad del 80% para la gestión interna de memoria de CD-HIT
MEM_MB=$(( RAW_MEM * 80 / 100 ))

# Hilos = CPUs reservadas (el lanzador ya pregunta el número, por defecto 8)
THREADS="${SLURM_CPUS_PER_TASK}"

# ------------------------------------------------------------------------------
# 5. EJECUCIÓN DE CD-HIT
# ------------------------------------------------------------------------------
echo "=========================================================="
echo " EJECUCIÓN DE CD-HIT EN NODO DE CÓMPUTO"
echo "=========================================================="
echo " -> Log en tiempo real: $RUN_LOG"
echo "    (síguelo con: tail -f $RUN_LOG)"
echo " -> Job ID Slurm:       ${SLURM_JOB_ID}"
echo " -> Nodo asignado:      $(hostname)"
echo " -> Herramienta:        $CDHIT_CMD"
echo " -> Hilos asignados:    $THREADS (de ${SLURM_CPUS_PER_TASK} reservados)"
echo " -> Límite RAM CD-HIT:  ${MEM_MB} MB"
echo " -> Archivo entrada:    $INPUT_FASTA"
echo " -> Carpeta resultados: $DEST_DIR"
echo " -> Banderas aplicadas: $CDHIT_EXTRA"
echo "=========================================================="

# En segundo plano + wait: así el trap USR1 se atiende mientras CD-HIT corre
set +e
"$CDHIT_CMD" -i input.fasta -o output_cdhit -T "$THREADS" -M "$MEM_MB" $CDHIT_EXTRA > "$RUN_LOG" 2>&1 &
TOOL_PID=$!
CDHIT_EXIT_CODE=0
while kill -0 "$TOOL_PID" 2>/dev/null; do
    wait "$TOOL_PID"; CDHIT_EXIT_CODE=$?
done
set -e

if [ $CDHIT_EXIT_CODE -ne 0 ]; then
    echo " [ERROR CRÍTICO] CD-HIT falló con código de salida $CDHIT_EXIT_CODE." >&2
    if [ -f "$RUN_LOG" ]; then
        echo " --- Últimas líneas del log de CD-HIT ($RUN_LOG) ---" >&2
        tail -n 30 "$RUN_LOG" >&2
    fi
    exit $CDHIT_EXIT_CODE
fi

if [ ! -s output_cdhit ]; then
    echo " [ERROR CRÍTICO] El archivo FASTA generado por CD-HIT está vacío." >&2
    exit 1
fi

# ------------------------------------------------------------------------------
# 6. PARSEO Y GENERACIÓN DEL REPORTE DE HAPLOTIPOS
# ------------------------------------------------------------------------------
echo " Generando informe estructurado de haplotipos..."

awk '
/^>Cluster/ {
    if (rep != "") {
        print "HAPLOTIPO REPRESENTATIVO: " rep
        print "Total secuencias colapsadas en este grupo: " count
        print "Miembros agrupados:"
        print (seqs != "" ? seqs : "  (Ninguno, secuencia única)")
        print "--------------------------------------------------------------------------------"
    }
    rep=""
    seqs=""
    count=0
    next
}
{
    match($0, />[^ ]+/)
    seq_id = substr($0, RSTART+1, RLENGTH-1)
    sub(/\.\.\.$/, "", seq_id)
    
    match($0, /at [0-9\.]+%/)
    ident = (RLENGTH > 0) ? substr($0, RSTART, RLENGTH) : "100.00% (referencia)"

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
        print "Miembros agrupados:"
        print (seqs != "" ? seqs : "  (Ninguno, secuencia única)")
        print "--------------------------------------------------------------------------------"
    }
}' output_cdhit.clstr > haplotypes_summary.txt

# ------------------------------------------------------------------------------
# 7. CONSOLIDACIÓN DE RESULTADOS
# ------------------------------------------------------------------------------
cp output_cdhit "$DEST_DIR/${BASENAME}_representatives.fasta"
cp output_cdhit.clstr "$DEST_DIR/${BASENAME}_raw.clstr"
cp haplotypes_summary.txt "$DEST_DIR/${BASENAME}_haplotypes.txt"
COMPLETADO=1

NUM_IN=$(grep -c '^>' input.fasta || true)
NUM_REPS=$(grep -c '^>' output_cdhit || true)

# Registro de metadatos de ejecución
{
    echo "=========================================================="
    echo " REGISTRO DE EJECUCIÓN - CD-HIT DEDUPLICATION"
    echo "=========================================================="
    echo " Fecha de ejecución:     $(date '+%Y-%m-%d %H:%M:%S')"
    echo " Slurm Job ID:           ${SLURM_JOB_ID}"
    echo " Archivo original:       $INPUT_FASTA"
    echo " Hash MD5 entrada:       $(md5sum "$INPUT_FASTA" | awk '{print $1}')"
    echo " Herramienta:            $CDHIT_CMD"
    echo " Secuencias entrada:     $NUM_IN"
    echo " Haplotipos/Clusters:    $NUM_REPS"
    echo " Hilos usados:           $THREADS"
    echo " Directorio de salida:   $DEST_DIR"
    echo " Log de ejecución:       $RUN_LOG"
    echo " Comando exacto:         $CDHIT_CMD -i input.fasta -o output_cdhit -T $THREADS -M $MEM_MB $CDHIT_EXTRA"
    echo "=========================================================="
} > "$INFO_FILE"

echo "=========================================================="
echo " Proceso completado con éxito."
echo " Resultados consolidados en: $DEST_DIR"
echo "   - FASTA de representativos: ${BASENAME}_representatives.fasta ($NUM_REPS haplotipos)"
echo "   - Informe de haplotipos:     ${BASENAME}_haplotypes.txt"
echo "   - Matriz cruda de clusters:  ${BASENAME}_raw.clstr"
echo "   - Registro de parámetros:    $INFO_FILE"
echo "=========================================================="