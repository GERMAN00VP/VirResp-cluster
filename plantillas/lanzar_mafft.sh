#!/bin/bash
# ==============================================================================
# LANZADOR INTERACTIVO DE MAFFT (VirResp-cluster)
# Uso: ./lanzar_mafft.sh alineamiento.fasta
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
echo " ALINEAMIENTO MULTIPLE DE SECUENCIAS CON MAFFT"
linea
echo " Fichero de entrada: Secuencias sin alinear FASTA (*.fasta, *.fa, *.fas)"
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
MASTER_SCRIPT="$(realpath -m "$SCRIPT_DIR/../scripts/master_mafft.sh")"

if [ ! -x "$MASTER_SCRIPT" ]; then
    error "No existe o no es ejecutable el script maestro:" "$MASTER_SCRIPT"
fi

# ------------------------------------------------------------------------------
# 2. INSPECCIÓN DEL ARCHIVO FASTA
# ------------------------------------------------------------------------------
read -r NUM_SEQS MIN_LEN MAX_LEN AVG_LEN < <(
    tr -d '\r' < "$INPUT_FASTA" | awk '
        BEGIN { n = 0; seen = 0; min = 0; max = 0; total = 0 }
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
        }
        END {
            if (n > 0) {
                if (!seen || len < min) { min = len; seen = 1 }
                if (len > max) max = len
                total += len
            }
            avg = (n > 0) ? int(total / n) : 0
            print n, min, max, avg
        }'
)

if [ "$NUM_SEQS" -eq 0 ]; then
    error "No se detectaron cabeceras FASTA (líneas '>') en '$INPUT_FASTA'."
fi

if [ "$NUM_SEQS" -lt 2 ]; then
    error "Se detectó solo $NUM_SEQS secuencia." "Se requieren al menos 2 secuencias para realizar un alineamiento."
fi

echo " INSPECCIÓN DE ENTRADA"
echo " -> Archivo:            $INPUT_FASTA"
echo " -> Secuencias totales: $NUM_SEQS"
echo " -> Rango de longitud:  $MIN_LEN a $MAX_LEN pb/aa (Promedio: ~$AVG_LEN)"

if [ "$MIN_LEN" -eq "$MAX_LEN" ] && [ "$MAX_LEN" -gt 0 ]; then
    echo ""
    echo " [AVISO] Todas las secuencias tienen exactamente la misma longitud ($MAX_LEN)."
    echo "         Este archivo podría estar ya alineado."
    if ! ask_yn "         ¿Deseas volver a alinear de todos modos?" "n"; then
        echo "Operación cancelada."
        exit 0
    fi
fi

DUPS=$(tr -d '\r' < "$INPUT_FASTA" | grep '^>' | sed 's/^>//; s/[[:space:]].*$//' | sort | uniq -d | head -n 5 || true)
if [ -n "$DUPS" ]; then
    echo ""
    echo " [AVISO] Se detectaron identificadores duplicados en el FASTA:"
    echo "         $(echo "$DUPS" | tr '\n' ' ')"
    echo "         MAFFT puede renombrar o confundir muestras con el mismo ID."
    if ! ask_yn "         ¿Continuar de todos modos?" "n"; then
        echo "Operación cancelada."
        exit 0
    fi
fi
linea

# ------------------------------------------------------------------------------
# 3. SELECCIÓN DE PARÁMETROS DE ALINEAMIENTO
# ------------------------------------------------------------------------------
echo ""
echo "1. ESTRATEGIA DE ALINEAMIENTO"
echo "  1) AUTO (FFT-NS-1, FFT-NS-2, FFT-NS-i o L-INS-i según el tamaño del conjunto)"
echo "  2) L-INS-i (Precisión local máxima; ideal para < 200 seqs con 1 dominio conservado)"
echo "  3) G-INS-i (Precisión global máxima; ideal para < 200 seqs con homología global)"
echo "  4) E-INS-i (Recomendado para < 200 seqs con múltiples dominios y grandes gaps)"
echo "  5) FFT-NS-2 (Rápido; método progresivo estándar para miles de secuencias)"
echo "  6) PartTree (Ultra-masivo; obligatorio para > 10.000 genomas o secuencias largas)"
ask_choice "Selecciona estrategia" "1" 6
STRAT_OPT="$CHOICE"

case "$STRAT_OPT" in
    1) STRAT_FLAGS="--auto";                              STRAT_TXT="AUTO (Automático según dataset)" ;;
    2) STRAT_FLAGS="--localpair --maxiterate 1000";       STRAT_TXT="L-INS-i (Precisión Local)" ;;
    3) STRAT_FLAGS="--globalpair --maxiterate 1000";      STRAT_TXT="G-INS-i (Precisión Global)" ;;
    4) STRAT_FLAGS="--genafpair --maxiterate 1000";       STRAT_TXT="E-INS-i (Especial para Gaps Complejos)" ;;
    5) STRAT_FLAGS="--retree 2";                          STRAT_TXT="FFT-NS-2 (Rápido Progresivo)" ;;
    6) STRAT_FLAGS="--parttree --retree 1";               STRAT_TXT="PartTree (Ultra-masivo)" ;;
esac

echo ""
echo "2. AJUSTE DE ORIENTACIÓN DE HEBRAS (--adjustdirection)"
echo "   Útil en genomas/fragmentos de ADN/ARN donde algunas secuencias están en hebra reversa."
if ask_yn "   ¿Detectar y reorientar automáticamente secuencias invertidas?" "s"; then
    ADJUST_FLAG="--adjustdirection"
    ADJUST_TXT="Activado (--adjustdirection)"
else
    ADJUST_FLAG=""
    ADJUST_TXT="Desactivado (conservar hebra original)"
fi

echo ""
echo "3. TRATAMIENTO DE REGIONES HIPERVARIABLES E INSERCIONES NO HOMÓLOGAS"
echo "   Explicación para genomas virales:"
echo "   Cuando los datos genómicos contienen inserciones únicas, deleciones grandes o extremos"
echo "   incompletos, los métodos globales intentan forzar un alineamiento espurio de esas zonas."
echo "   La opción 'Leave gappy regions' deja estas regiones no alineadas (en bloques de gaps),"
echo "   evitando distorsionar los bloques homólogos bien conservados."
echo ""
echo "  1) Desactivado (Alineación uniforme continua predeterminada)"
echo "  2) Activado: Leave gappy regions con Unalignlevel = 0.8 (Recomendado para genomas virales)"
echo "  3) Personalizado (Ajustar Unalignlevel entre 0.1 y 0.8)"
ask_choice "Selección" "1" 3
GAPPY_OPT="$CHOICE"

case "$GAPPY_OPT" in
    1)
        GAPPY_FLAGS=""
        GAPPY_TXT="Desactivado (alineamiento global estándar)"
        ;;
    2)
        GAPPY_FLAGS="--leavegappyregion --unalignlevel 0.8"
        GAPPY_TXT="Activado (--leavegappyregion --unalignlevel 0.8)"
        ;;
    3)
        while true; do
            read -r -p "   Introduce el valor de Unalignlevel (0.1 a 0.8) [Intro = 0.8]: " UNALIGN_VAL
            UNALIGN_VAL="${UNALIGN_VAL:-0.8}"
            if [[ "$UNALIGN_VAL" =~ ^0\.[1-8]$ ]]; then
                break
            fi
            echo "   Valor no válido. Debe ser un decimal entre 0.1 y 0.8 (ej. 0.5)."
        done
        GAPPY_FLAGS="--leavegappyregion --unalignlevel $UNALIGN_VAL"
        GAPPY_TXT="Activado (--leavegappyregion --unalignlevel $UNALIGN_VAL)"
        ;;
esac

if [ -n "$GAPPY_FLAGS" ]; then
    if [ "$NUM_SEQS" -gt 1000 ]; then
        echo ""
        echo " [AVISO] MAFFT desaconseja 'Leave gappy regions' para conjuntos > 1.000 secuencias"
        echo "         debido a la alta complejidad de cómputo O(N^2)."
        if ! ask_yn "         ¿Deseas mantener esta opción de todos modos?" "n"; then
            GAPPY_FLAGS=""
            GAPPY_TXT="Desactivado (por límite de secuencias)"
        fi
    fi
    # Si se selecciona leavegappyregion y la estrategia es progresiva o PartTree, forzar cálculo de pares homólogos si no estaba fijado
    if [ "$STRAT_OPT" -eq 1 ] || [ "$STRAT_OPT" -eq 5 ]; then
        STRAT_FLAGS="--globalpair $STRAT_FLAGS"
        STRAT_TXT="$STRAT_TXT + Global Pairwise (requerido por --leavegappyregion)"
    fi
fi

echo ""
echo "4. REORDENACIÓN DE SECUENCIAS SEGÚN SIMILITUD (--reorder)"
if ask_yn "   ¿Reordenar las secuencias alineadas por similitud filogenética?" "n"; then
    REORDER_FLAG="--reorder"
    REORDER_TXT="Activado (--reorder)"
else
    REORDER_FLAG=""
    REORDER_TXT="Desactivado (mantener orden del FASTA original)"
fi

echo ""
echo "5. PENALIZACIÓN DE APERTURA DE GAPS (--op)"
read -r -p "   Valor de penalización de apertura de gaps [Intro = 1.53 estándar]: " OP_VAL
OP_VAL="${OP_VAL:-1.53}"
if [ "$OP_VAL" != "1.53" ]; then
    OP_FLAG="--op $OP_VAL"
else
    OP_FLAG=""
fi

# ------------------------------------------------------------------------------
# 4. CONSTRUCCIÓN DE BANDERAS Y ESTIMACIÓN DE RECURSOS SLURM
# ------------------------------------------------------------------------------
FLAGS=()
FLAGS+=($STRAT_FLAGS)

if [ -n "$ADJUST_FLAG" ]; then
    FLAGS+=("$ADJUST_FLAG")
fi

if [ -n "$GAPPY_FLAGS" ]; then
    FLAGS+=($GAPPY_FLAGS)
fi

if [ -n "$REORDER_FLAG" ]; then
    FLAGS+=("$REORDER_FLAG")
fi

if [ -n "$OP_FLAG" ]; then
    FLAGS+=($OP_FLAG)
fi

MAFFT_EXTRA="${FLAGS[*]}"

# Cálculo de recursos Slurm
if [ "$STRAT_OPT" -eq 6 ] || [ "$NUM_SEQS" -gt 10000 ]; then
    PARTITION="long_idx";   TIME="5-00:00:00"; CPUS=32; MEM="128G"
elif [ "$STRAT_OPT" -in_set "2 3 4" ] || [ -n "$GAPPY_FLAGS" ] || [ "$NUM_SEQS" -gt 1000 ]; then
    PARTITION="middle_idx"; TIME="48:00:00";   CPUS=16; MEM="64G"
else
    PARTITION="short_idx";  TIME="12:00:00";   CPUS=16; MEM="32G"
fi

# ------------------------------------------------------------------------------
# 5. RESUMEN Y CONFIRMACIÓN
# ------------------------------------------------------------------------------
echo ""
linea
echo " RESUMEN DE LA CONFIGURACIÓN DE MAFFT"
linea
echo " Archivo de entrada:     $INPUT_FASTA ($NUM_SEQS secuencias)"
echo " Estrategia principal:   $STRAT_TXT"
echo " Orientación de hebras:  $ADJUST_TXT"
echo " Regiones gappy/hiper:   $GAPPY_TXT"
echo " Orden de salida:        $REORDER_TXT"
echo " Banderas finales:       mafft --thread \$SLURM_CPUS_PER_TASK $MAFFT_EXTRA"
echo " Salida esperada:        $WORKDIR/$(basename "$INPUT_FASTA" | sed -E 's/\.(fasta|fa|fas)$//')_aligned.fasta"
echo ""
echo " Reserva en Slurm:"
echo "   - Cómputo:  $CPUS CPUs | Memoria: $MEM"
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
    --output="$WORKDIR/mafft_%j.out" \
    --error="$WORKDIR/mafft_%j.err" \
    --partition="$PARTITION" \
    --time="$TIME" \
    --cpus-per-task="$CPUS" \
    --mem="$MEM" \
    --job-name="MAFFT_${NUM_SEQS}" \
    "$MASTER_SCRIPT" "$INPUT_FASTA" "$MAFFT_EXTRA")
JOB_ID="${JOB_ID%%;*}"

echo ""
echo " Job enviado con éxito a Slurm. ID del trabajo: $JOB_ID"
echo " Seguimiento y gestión:"
echo "   - Estado de la cola:            squeue -u \$USER"
echo "   - Ver salida en tiempo real:     tail -f $WORKDIR/mafft_${JOB_ID}.err"