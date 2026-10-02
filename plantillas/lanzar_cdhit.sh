#!/bin/bash
# ==============================================================================
# LANZADOR INTERACTIVO DE CD-HIT (VirResp-cluster)
# Uso: ./lanzar_cdhit.sh mi_fichero.fasta
# ==============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# Funciones auxiliares
# ------------------------------------------------------------------------------
linea() { echo "=========================================================="; }

error() {
    echo ""
    echo " [ERROR] $1"
    shift
    for extra in "$@"; do echo "         $extra"; done
    echo ""
    exit 1
}

ask_choice() {
    local ans
    while true; do
        read -r -p "$1 [1-$3] (Intro = $2): " ans || { echo ""; echo "Entrada cancelada."; exit 1; }
        ans="${ans:-$2}"
        if [[ "$ans" =~ ^[0-9]+$ ]] && [ "$ans" -ge 1 ] && [ "$ans" -le "$3" ]; then
            CHOICE="$ans"
            return 0
        fi
        echo "   Opción no válida: escribe un número entre 1 y $3."
    done
}

ask_yn() {
    local ans def="$2" hint
    if [ "$def" = "s" ]; then hint="S/n"; else hint="s/N"; fi
    while true; do
        read -r -p "$1 [$hint]: " ans || { echo ""; echo "Entrada cancelada."; exit 1; }
        ans="${ans:-$def}"
        case "$ans" in
            s|S|si|SI|Si|sí|Sí|y|Y) return 0 ;;
            n|N|no|NO|No)           return 1 ;;
        esac
        echo "   Responde 's' (sí) o 'n' (no)."
    done
}

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

if ! command -v sbatch >/dev/null 2>&1; then
    error "No se encuentra el gestor Slurm (sbatch)." "Ejecuta este script dentro del clúster."
fi

INPUT_FASTA="$(realpath "$1")"
WORKDIR="$(dirname "$INPUT_FASTA")"

if [ ! -w "$WORKDIR" ]; then
    error "Sin permiso de escritura en la carpeta:" "$WORKDIR"
fi

SCRIPT_DIR="$(cd "$(dirname "$(realpath "${BASH_SOURCE[0]}")")" && pwd)"
MASTER_SCRIPT="$(realpath -m "$SCRIPT_DIR/../scripts/master_cdhit.sh")"

if [ ! -x "$MASTER_SCRIPT" ]; then
    error "No existe o no es ejecutable el script maestro:" "$MASTER_SCRIPT"
fi

# Sanitizar saltos de línea (\r)
sed -i 's/\r$//' "$INPUT_FASTA"

# ------------------------------------------------------------------------------
# 2. INSPECCIÓN DEL ARCHIVO FASTA
# ------------------------------------------------------------------------------
read -r NUM_SEQS MIN_LEN MAX_LEN AVG_LEN TOTAL_CHARS NUC_CHARS < <(
    tr -d '\r' < "$INPUT_FASTA" | awk '
        BEGIN { n = 0; seen = 0; min = 0; max = 0; total = 0; nuc = 0 }
        /^>/ {
            if (n > 0) {
                if (!seen || len < min) { min = len; seen = 1 }
                if (len > max) max = len
                total += len
            }
            n++; len = 0; next
        }
        {
            gsub(/[ \t]/, "")
            len += length($0)
            s = toupper($0)
            t = s; nuc += gsub(/[ACGTUN]/, "", t)
        }
        END {
            if (n > 0) {
                if (!seen || len < min) { min = len; seen = 1 }
                if (len > max) max = len
                total += len
            }
            avg = (n > 0) ? int(total / n) : 0
            print n, min, max, avg, total + 0, nuc + 0
        }'
)

if [ "$NUM_SEQS" -eq 0 ]; then
    error "No se detectaron cabeceras FASTA (líneas '>') en '$INPUT_FASTA'."
fi

echo " INSPECCIÓN DE ENTRADA"
echo " -> Archivo:            $INPUT_FASTA"
echo " -> Secuencias totales: $NUM_SEQS"
echo " -> Longitud estimada:  $MIN_LEN a $MAX_LEN pb/aa (Promedio: ~$AVG_LEN)"
linea

# Autodetección del tipo de secuencia (ADN/ARN vs Proteína)
PCT_NUC=0
if [ "$TOTAL_CHARS" -gt 0 ]; then
    PCT_NUC=$(( 100 * NUC_CHARS / TOTAL_CHARS ))
fi

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
echo "   Autodetección basada en composición: $DETECTED_LABEL"
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
CPUS=8
MEM="32G"

if [ "$NUM_SEQS" -gt 50000 ]; then
    PARTITION="middle_idx"
    TIME="24:00:00"
    CPUS=16
    MEM="64G"
fi

BASENAME="$(basename "$INPUT_FASTA")"
BASENAME="${BASENAME%.fasta}"; BASENAME="${BASENAME%.fa}"
BASENAME="${BASENAME%.fas}";   BASENAME="${BASENAME%.aligned}"

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
echo " Comand line flags:      $CDHIT_CMD -i input.fasta -o output_cdhit $CDHIT_EXTRA"
echo " Directorio salida:      $WORKDIR/${BASENAME}_haplotypes/"
echo ""
echo " Reserva en Slurm:"
echo "   - Cómputo:  $CPUS CPUs | Memoria RAM: $MEM"
echo "   - Tiempo:   máximo $TIME | Cola: $PARTITION"
linea

if ! ask_yn "¿Confirmar y enviar trabajo a Slurm?" "s"; then
    echo "Operación cancelada por el usuario."
    exit 0
fi

# ------------------------------------------------------------------------------
# 6. ENVÍO DEL JOB A SLURM
# ------------------------------------------------------------------------------
JOB_ID=$(sbatch --parsable \
    --chdir="$WORKDIR" \
    --output="$WORKDIR/cdhit_%j.out" \
    --error="$WORKDIR/cdhit_%j.err" \
    --partition="$PARTITION" \
    --time="$TIME" \
    --cpus-per-task="$CPUS" \
    --mem="$MEM" \
    --job-name="CDHIT_${NUM_SEQS}" \
    "$MASTER_SCRIPT" "$INPUT_FASTA" "$CDHIT_CMD" "$CDHIT_EXTRA")
JOB_ID="${JOB_ID%%;*}"

echo ""
echo " Job enviado con éxito a Slurm. ID del trabajo: $JOB_ID"
echo ""
echo " Seguimiento y gestión:"
echo "   - Estado de la cola:            squeue -u \$USER"
echo "   - Ver salida en tiempo real:     tail -f $WORKDIR/cdhit_${JOB_ID}.err"
echo "   - Directorio de resultados:     $WORKDIR/${BASENAME}_haplotypes/"