#!/bin/bash
# ==============================================================================
# LANZADOR INTERACTIVO DE CD-HIT (VirResp-cluster)
# Uso: ./lanzar_cdhit.sh mi_fichero.fasta
# ==============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# Funciones comunes (linea, error, ask_choice, ask_yn, comprobaciones de FASTA)
# ------------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "$(realpath "${BASH_SOURCE[0]}")")" && pwd)"
# shellcheck source=comun.sh
source "$SCRIPT_DIR/comun.sh"

# ------------------------------------------------------------------------------
# Encabezado
# ------------------------------------------------------------------------------
linea
echo " DESDUPLICACIÓN Y AGRUPAMIENTO DE SECUENCIAS CON CD-HIT"
linea
echo " Fichero de entrada: Secuencias nucleotídicas o aminoácidos en FASTA"
linea

# ------------------------------------------------------------------------------
# 1. COMPROBACIONES PREVIAS
# ------------------------------------------------------------------------------
if [ "$#" -ne 1 ]; then
    error "Indica el archivo FASTA de entrada." "Uso: $0 <secuencias.fasta>"
fi

if [ ! -f "$1" ]; then
    error "El archivo '$1' no existe."
fi

if [ ! -s "$1" ]; then
    error "El archivo '$1' está vacío."
fi

case "$1" in
    *.gz) error "El archivo está comprimido (.gz)." "Descomprímelo con: gunzip $1" ;;
esac

comprobar_sbatch

INPUT_FASTA="$(realpath "$1")"
WORKDIR="$(dirname "$INPUT_FASTA")"

if [ ! -w "$WORKDIR" ]; then
    error "Sin permiso de escritura en la carpeta:" "$WORKDIR"
fi

MASTER_SCRIPT="$(realpath -m "$SCRIPT_DIR/../scripts/master_cdhit.sh")"

if [ ! -x "$MASTER_SCRIPT" ]; then
    error "No existe o no es ejecutable el script maestro:" "$MASTER_SCRIPT"
fi

# ------------------------------------------------------------------------------
# 2. INSPECCIÓN DEL ARCHIVO FASTA
# ------------------------------------------------------------------------------
fasta_basico "$INPUT_FASTA"

if [ "$NUM_SEQS" -eq 0 ]; then
    error "No se detectaron cabeceras FASTA (líneas '>') en '$INPUT_FASTA'."
fi

fasta_longitudes_muestra "$INPUT_FASTA" 50

echo " INSPECCIÓN DE ENTRADA"
echo " -> Archivo:            $INPUT_FASTA"
echo " -> Secuencias totales: $NUM_SEQS"
echo " -> Longitud estimada:  $MIN_LEN a $MAX_LEN pb/aa (Promedio: ~$AVG_LEN; muestra de $N_MUESTRA primeras)"
linea

# Autodetección del tipo de secuencia (ADN/ARN vs Proteína) sobre una muestra
fasta_tipo_muestra "$INPUT_FASTA"
PCT_NUC="$NT_PCT"

if [ "$PCT_NUC" -ge 85 ]; then
    DEFAULT_SEQ_TYPE=1
    DETECTED_LABEL="Nucleótidos (ADN/ARN) -> cd-hit-est"
else
    DEFAULT_SEQ_TYPE=2
    DETECTED_LABEL="Aminoácidos (Proteínas) -> cd-hit"
fi

# ------------------------------------------------------------------------------
# 3. SELECCIÓN DE PARÁMETROS
# ------------------------------------------------------------------------------
echo ""
echo "1. TIPO DE SECUENCIA"
echo "   Autodetección basada en composición (muestra de ~2 MB): $DETECTED_LABEL"
echo "  1) Nucleótidos / ADN / ARN (cd-hit-est)"
echo "  2) Aminoácidos / Proteínas (cd-hit)"
ask_choice "Selección" "$DEFAULT_SEQ_TYPE" 2

if [ "$CHOICE" -eq 1 ]; then
    CDHIT_CMD="cd-hit-est"
    SEQ_TYPE_TXT="Nucleótidos (cd-hit-est)"
else
    CDHIT_CMD="cd-hit"
    SEQ_TYPE_TXT="Proteínas (cd-hit)"
fi

echo ""
echo "2. UMBRAL DE SIMILITUD DE IDENTIDAD (-c)"
echo "  1) 100% Identidad (-c 1.0)  [Desduplicación pura / Haplotipado estricto]"
echo "  2)  99% Identidad (-c 0.99) [Colapsar errores de PCR/secuenciación o microvariantes]"
echo "  3)  95% Identidad (-c 0.95) [Agrupar por linajes/subclados virales]"
echo "  4) Personalizado            [Introducir valor numérico]"
ask_choice "Selección" "1" 4

case "$CHOICE" in
    1) ID_THRESHOLD="1.0" ;;
    2) ID_THRESHOLD="0.99" ;;
    3) ID_THRESHOLD="0.95" ;;
    4)
        while true; do
            read -r -p "   Introduce el umbral de identidad (0.80 a 1.0) [Intro = 1.0]: " ID_INPUT
            ID_INPUT="${ID_INPUT:-1.0}"
            if awk "BEGIN{exit !($ID_INPUT >= 0.80 && $ID_INPUT <= 1.0)}"; then
                ID_THRESHOLD="$ID_INPUT"
                break
            fi
            echo "   Valor no válido. Debe estar entre 0.80 y 1.0 (ej. 0.92)."
        done
        ;;
esac

# Ajuste automático e infalible del tamaño de palabra (-n) según requerimientos internos de CD-HIT
if [ "$CDHIT_CMD" = "cd-hit-est" ]; then
    if awk "BEGIN{exit !($ID_THRESHOLD >= 0.95)}"; then WORD_SIZE=10
    elif awk "BEGIN{exit !($ID_THRESHOLD >= 0.90)}"; then WORD_SIZE=8
    elif awk "BEGIN{exit !($ID_THRESHOLD >= 0.88)}"; then WORD_SIZE=6
    elif awk "BEGIN{exit !($ID_THRESHOLD >= 0.85)}"; then WORD_SIZE=5
    else WORD_SIZE=4; fi
else
    if awk "BEGIN{exit !($ID_THRESHOLD >= 0.70)}"; then WORD_SIZE=5
    elif awk "BEGIN{exit !($ID_THRESHOLD >= 0.60)}"; then WORD_SIZE=4
    elif awk "BEGIN{exit !($ID_THRESHOLD >= 0.50)}"; then WORD_SIZE=3
    else WORD_SIZE=2; fi
fi

echo ""
echo "3. COBERTURA MÍNIMA DE ALINEAMIENTO (-aS)"
echo "   Define qué porcentaje de la secuencia más corta debe estar cubierta en el agrupamiento."
read -r -p "   Valor de cobertura (0.0 a 1.0) [Intro = 1.0 (100%)]: " AS_VAL
AS_VAL="${AS_VAL:-1.0}"

# Parámetros fijos recomendados:
# -g 1: agrupar en el cluster más similar
# -d 0: no truncar la cabecera FASTA en el reporte .clstr
CDHIT_EXTRA="-c $ID_THRESHOLD -n $WORD_SIZE -aS $AS_VAL -g 1 -d 0"

# ------------------------------------------------------------------------------
# 4. ESTIMACIÓN DE RECURSOS SLURM
# ------------------------------------------------------------------------------
PARTITION="short_idx"
TIME="04:00:00"
MEM="32G"

if [ "$NUM_SEQS" -gt 50000 ]; then
    PARTITION="middle_idx"
    TIME="24:00:00"
    MEM="64G"
fi

pedir_hilos "4. HILOS DE CÓMPUTO (-T)" 8 \
    "CD-HIT apenas escala por encima de 8 hilos." \
    "Se reservarán tantas CPUs como hilos."
CPUS="$HILOS"

BASENAME="$(basename "$INPUT_FASTA")"
BASENAME="${BASENAME%.fasta}"; BASENAME="${BASENAME%.fa}"
BASENAME="${BASENAME%.fas}";   BASENAME="${BASENAME%.aligned}"

pedir_outdir "$WORKDIR/${BASENAME}_haplotypes" "$WORKDIR"

# ------------------------------------------------------------------------------
# 5. RESUMEN Y CONFIRMACIÓN
# ------------------------------------------------------------------------------
echo ""
linea
echo " RESUMEN DE LA CONFIGURACIÓN DE CD-HIT"
linea
echo " Archivo de entrada:     $INPUT_FASTA ($NUM_SEQS secuencias)"
echo " Herramienta / Modo:     $SEQ_TYPE_TXT"
echo " Umbral de identidad:    $ID_THRESHOLD (Word size: $WORD_SIZE)"
echo " Cobertura mínima:       $AS_VAL"
echo " Comand line flags:      $CDHIT_CMD -i input.fasta -o output_cdhit -T $CPUS $CDHIT_EXTRA"
echo " Directorio salida:      $OUTDIR/"
echo ""
echo " Reserva en Slurm:"
echo "   - Cómputo:  $CPUS CPUs (= hilos de CD-HIT) | Memoria RAM: $MEM"
echo "   - Tiempo:   máximo $TIME | Cola: $PARTITION"
linea

if ! ask_yn "¿Confirmar y enviar trabajo a Slurm?" "s"; then
    echo "Operación cancelada por el usuario."
    exit 0
fi

# ------------------------------------------------------------------------------
# 6. ENVÍO DEL JOB A SLURM
# ------------------------------------------------------------------------------
enviar_job \
    --chdir="$WORKDIR" \
    --output="$OUTDIR/slurm_%j.out" \
    --error="$OUTDIR/slurm_%j.err" \
    --partition="$PARTITION" \
    --time="$TIME" \
    --cpus-per-task="$CPUS" \
    --mem="$MEM" \
    --signal=B:USR1@1800 \
    --job-name="CDHIT_${NUM_SEQS}" \
    "$MASTER_SCRIPT" "$INPUT_FASTA" "$CDHIT_CMD" "$CDHIT_EXTRA" "$OUTDIR"

echo ""
echo " Job enviado con éxito a Slurm. ID del trabajo: $JOB_ID"
echo ""
echo " Seguimiento y gestión:"
echo "   - Estado de la cola:            squeue -u \$USER"
echo "   - Ver salida en tiempo real:     tail -f $OUTDIR/slurm_${JOB_ID}.err"
echo "   - Directorio de resultados:     $OUTDIR/"