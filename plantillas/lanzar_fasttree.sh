#!/bin/bash
# ==============================================================================
# LANZADOR INTERACTIVO DE FASTTREE (VirResp-cluster)
# Uso: ./lanzar_fasttree.sh alineamiento.fasta
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
echo " ÁRBOL FILOGENÉTICO DE MÁXIMA VEROSIMILITUD CON FASTTREE"
linea
echo " Fichero de entrada: Alineamiento FASTA (*.fasta, *.aln, *.fa)"
linea

# ------------------------------------------------------------------------------
# 1. COMPROBACIONES PREVIAS
# ------------------------------------------------------------------------------
if [ "$#" -ne 1 ]; then
    error "Indica el alineamiento de entrada." "Uso: $0 <alineamiento.fasta>"
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
MASTER_SCRIPT="$(realpath -m "$SCRIPT_DIR/../scripts/master_fasttree.sh")"

if [ ! -x "$MASTER_SCRIPT" ]; then
    error "No existe o no es ejecutable el script maestro:" "$MASTER_SCRIPT"
fi

# ------------------------------------------------------------------------------
# 2. INSPECCIÓN DEL ALINEAMIENTO
# ------------------------------------------------------------------------------
read -r NUM_SEQS MIN_LEN MAX_LEN TOTAL_CHARS NUC_CHARS GAP_CHARS < <(
    tr -d '\r' < "$INPUT_FASTA" | awk '
        BEGIN { n = 0; seen = 0; min = 0; max = 0 }
        /^>/ {
            if (n > 0) {
                if (!seen || len < min) { min = len; seen = 1 }
                if (len > max) max = len
            }
            n++; len = 0; next
        }
        {
            gsub(/[ \t]/, "")
            len += length($0)
            s = toupper($0); total += length(s)
            t = s; nuc += gsub(/[ACGTUN]/, "", t)
            t = s; gap += gsub(/[-.?]/, "", t)
        }
        END {
            if (n > 0) {
                if (!seen || len < min) { min = len; seen = 1 }
                if (len > max) max = len
            }
            print n, min, max, total + 0, nuc + 0, gap + 0
        }'
)

if [ "$NUM_SEQS" -eq 0 ]; then
    error "No se detectaron cabeceras FASTA (líneas '>') en '$INPUT_FASTA'."
fi

if [ "$NUM_SEQS" -lt 4 ]; then
    error "Se detectaron $NUM_SEQS secuencias." "Se requieren al menos 4 secuencias para construir un árbol."
fi

ALN_LEN="$MAX_LEN"

echo " INSPECCIÓN DEL ALINEAMIENTO"
echo " -> Archivo:      $INPUT_FASTA"
echo " -> Secuencias:   $NUM_SEQS"
echo " -> Longitud:     $ALN_LEN pb/aa"

if [ "$MIN_LEN" -ne "$MAX_LEN" ]; then
    error "Longitudes desiguales detectadas ($MIN_LEN a $MAX_LEN posiciones)." \
          "El archivo no está alineado. Ejecuta primero 'lanzar_mafft.sh'."
fi

DUPS=$(tr -d '\r' < "$INPUT_FASTA" | grep '^>' | sed 's/^>//; s/[[:space:]].*$//' | sort | uniq -d | head -n 5 || true)
if [ -n "$DUPS" ]; then
    error "Identificadores de secuencia duplicados detectados:" "$(echo "$DUPS" | tr '\n' ' ')"
fi

BAD_IDS=$(tr -d '\r' < "$INPUT_FASTA" | grep -c -E '^>[^[:space:]]*[][(),:;]' || true)
if [ "$BAD_IDS" -gt 0 ]; then
    echo ""
    echo " [AVISO] $BAD_IDS nombres contienen caracteres no recomendados en Newick:  ( ) , : ; [ ]"
    if ! ask_yn "         ¿Deseas continuar de todos modos?" "n"; then
        echo "Operación cancelada."
        exit 0
    fi
fi
linea

# ------------------------------------------------------------------------------
# 3. SELECCIÓN DE PARÁMETROS FILOGENÉTICOS
# ------------------------------------------------------------------------------
RESIDUES=$(( TOTAL_CHARS - GAP_CHARS ))
if [ "$RESIDUES" -le 0 ]; then
    error "El alineamiento solo contiene huecos o caracteres ambiguos (-)."
fi
NT_PCT=$(( 100 * NUC_CHARS / RESIDUES ))

if [ "$NT_PCT" -ge 90 ]; then
    DETECTED=1
    DETECTED_TXT="Nucleótidos (ADN/ARN)"
else
    DETECTED=2
    DETECTED_TXT="Aminoácidos (Proteínas)"
fi

echo ""
echo "1. TIPO DE SECUENCIA"
echo "   Autodetección basada en composición: $DETECTED_TXT"
echo "  1) Nucleótidos (ADN / ARN)"
echo "  2) Aminoácidos (Proteínas)"
ask_choice "Selección" "$DETECTED" 2
SEQ_OPT="$CHOICE"

if [ "$SEQ_OPT" -eq 1 ]; then
    SEQ_TYPE="NT"
    SEQ_TYPE_TXT="Nucleótidos (ADN/ARN)"
    
    echo ""
    echo "2. MODELO DE SUSTITUCIÓN NUCLEOTÍDICA"
    echo "  1) GTR+CAT (General Time Reversible) [Recomendado para genomas/genes virales]"
    echo "  2) JC+CAT (Jukes-Cantor)"
    ask_choice "Selección" "1" 2
    NT_MODEL_OPT="$CHOICE"
    if [ "$NT_MODEL_OPT" -eq 1 ]; then
        MODEL_FLAG="-gtr"
        MODEL_TXT="GTR (General Time Reversible)"
    else
        MODEL_FLAG="" # JC es el valor por defecto de FastTree si no se pasa -gtr
        MODEL_TXT="JC (Jukes-Cantor)"
    fi
else
    SEQ_TYPE="AA"
    SEQ_TYPE_TXT="Aminoácidos (Proteínas)"
    
    echo ""
    echo "2. MODELO DE SUSTITUCIÓN EN AMINOÁCIDOS"
    echo "  1) LG (Le y Gascuel) [Recomendado para la mayoría de proteínas virales]"
    echo "  2) WAG (Whelan y Goldman)"
    echo "  3) JTT (Jones-Taylor-Thornton)"
    ask_choice "Selección" "1" 3
    AA_MODEL_OPT="$CHOICE"
    case "$AA_MODEL_OPT" in
        1) MODEL_FLAG="-lg";  MODEL_TXT="LG (Le-Gascuel)" ;;
        2) MODEL_FLAG="-wag"; MODEL_TXT="WAG (Whelan-Goldman)" ;;
        3) MODEL_FLAG="";     MODEL_TXT="JTT (Jones-Taylor-Thornton)" ;;
    esac
fi

echo ""
echo "3. NÚMERO DE CATEGORÍAS DE TASA DE EVOLUCIÓN (-cat)"
echo "   FastTree aproxima la variación de tasa entre sitios con el modelo CAT."
echo "  1) 20 categorías [Estándar por defecto en FastTree]"
echo "  2) 12 categorías [Ahorro de memoria/tiempo en datasets masivos]"
echo "  3) 8 categorías"
echo "  4) Personalizado"
ask_choice "Selección" "1" 4
CAT_OPT="$CHOICE"

case "$CAT_OPT" in
    1) CAT_NUM=20 ;;
    2) CAT_NUM=12 ;;
    3) CAT_NUM=8 ;;
    4)
        while true; do
            read -r -p "   Introduce el número de categorías CAT (4-50) [Intro = 20]: " CUSTOM_CAT
            CUSTOM_CAT="${CUSTOM_CAT:-20}"
            if [[ "$CUSTOM_CAT" =~ ^[0-9]+$ ]] && [ "$CUSTOM_CAT" -ge 4 ] && [ "$CUSTOM_CAT" -le 50 ]; then
                CAT_NUM="$CUSTOM_CAT"
                break
            fi
            echo "   Valor no válido. Debe ser un entero entre 4 y 50."
        done
        ;;
esac
CAT_FLAG="-cat $CAT_NUM"

echo ""
echo "4. OPTIMIZACIÓN LIKELIHOOD BAJO MODELO GAMMA (-gamma)"
echo "   Reoptimiza la verosimilitud (lnL) y las longitudes de rama bajo una distribución Gamma."
echo "  1) Sí (-gamma) [Recomendado para estimación precisa de longitudes de rama]"
echo "  2) No (mantener únicamente la aproximación CAT)"
ask_choice "Selección" "1" 2
GAMMA_OPT="$CHOICE"

if [ "$GAMMA_OPT" -eq 1 ]; then
    GAMMA_FLAG="-gamma"
    GAMMA_TXT="Activado (-gamma)"
else
    GAMMA_FLAG=""
    GAMMA_TXT="Desactivado (solo CAT)"
fi

echo ""
echo "5. ALGORITMO Y EXHAUSTIVIDAD DE BÚSQUEDA TOPOLÓGICA"
echo "  1) Búsqueda Estándar: NNI y SPR rápidos (equilibrio óptimo velocidad/precisión)"
echo "  2) Búsqueda Exhaustiva (-spr 4 -mlacc 2 -slownni): mayor profundidad en intercambios"
echo "     SPR y NNI con optimización ML más rigurosa (aumenta el tiempo x2 - x4)"
ask_choice "Selección" "1" 2
SEARCH_OPT="$CHOICE"

if [ "$SEARCH_OPT" -eq 2 ]; then
    SEARCH_FLAGS="-spr 4 -mlacc 2 -slownni"
    SEARCH_TXT="Exhaustiva (-spr 4 -mlacc 2 -slownni)"
else
    SEARCH_FLAGS=""
    SEARCH_TXT="Estándar (NNI/SPR por defecto)"
fi

echo ""
echo "6. SOPORTE DE RAMAS Y MÉTODOS DE EVALUACIÓN"
if [ "$NUM_SEQS" -le 5000 ]; then
    echo "  1) SH-like Local Supports (Shimodaira-Hasegawa) [Rápido, calculado por defecto]"
    echo "  2) Bootstrap Clásico (1000 repeticiones con -boot 1000) [Requiere mucho más tiempo]"
    echo "  3) Sin soporte de ramas (-nosupport)"
    ask_choice "Selección" "1" 3
    SUPPORT_OPT="$CHOICE"
else
    echo "  1) SH-like Local Supports (Shimodaira-Hasegawa) [Rápido, por defecto]"
    echo "  2) Sin soporte de ramas (-nosupport)"
    echo "   (Nota: Bootstrap con 1000 repeticiones desactivado para > 5.000 secuencias por coste de cómputo)"
    ask_choice "Selección" "1" 2
    if [ "$CHOICE" -eq 2 ]; then SUPPORT_OPT=3; else SUPPORT_OPT=1; fi
fi

case "$SUPPORT_OPT" in
    1) SUPPORT_FLAG="";          SUPPORT_TXT="SH-like Local Supports (Shimodaira-Hasegawa)" ;;
    2) SUPPORT_FLAG="-boot 1000"; SUPPORT_TXT="Bootstrap Clásico (1000 repeticiones)" ;;
    3) SUPPORT_FLAG="-nosupport"; SUPPORT_TXT="Sin soporte de ramas (-nosupport)" ;;
esac

# ------------------------------------------------------------------------------
# 4. CONSTRUCCIÓN DE COMANDOS Y ESTIMACIÓN DE RECURSOS
# ------------------------------------------------------------------------------
FLAGS=()

if [ "$SEQ_TYPE" = "NT" ]; then
    FLAGS+=(-nt)
fi

if [ -n "$MODEL_FLAG" ]; then
    FLAGS+=("$MODEL_FLAG")
fi

FLAGS+=($CAT_FLAG)

if [ -n "$GAMMA_FLAG" ]; then
    FLAGS+=("$GAMMA_FLAG")
fi

if [ -n "$SEARCH_FLAGS" ]; then
    FLAGS+=($SEARCH_FLAGS)
fi

if [ -n "$SUPPORT_FLAG" ]; then
    FLAGS+=($SUPPORT_FLAG)
fi

FLAGS+=(-seed 1253)

HUGE_NOTE=""
if [ "$NUM_SEQS" -gt 50000 ]; then
    FLAGS+=(-fastest)
    HUGE_NOTE="Aviso: Dataset masivo (>50.000 secuencias). Se añade automáticamente '-fastest' para optimización de memoria."
fi

FASTTREE_EXTRA="${FLAGS[*]}"

# Factor de esfuerzo para Slurm
WORK=$(( NUM_SEQS * ALN_LEN ))
EFF_WORK="$WORK"
if [ "$SEARCH_OPT" -eq 2 ];  then EFF_WORK=$(( EFF_WORK * 3 )); fi
if [ "$SUPPORT_OPT" -eq 2 ]; then EFF_WORK=$(( EFF_WORK * 10 )); fi

if [ "$NUM_SEQS" -gt 50000 ] || [ "$EFF_WORK" -gt 3000000000 ]; then
    PARTITION="long_idx";   TIME="5-00:00:00"; CPUS=16
elif [ "$NUM_SEQS" -gt 5000 ] || [ "$EFF_WORK" -gt 200000000 ]; then
    PARTITION="middle_idx"; TIME="48:00:00";   CPUS=8
else
    PARTITION="short_idx";  TIME="12:00:00";   CPUS=4
fi

if [ "$SEQ_TYPE" = "NT" ]; then BYTES_POS=16; else BYTES_POS=80; fi
EST_BYTES=$(( 3 * NUM_SEQS * ALN_LEN * BYTES_POS ))
EST_GB=$(( EST_BYTES / 1000000000 + 5 ))
EST_GB=$(( (EST_GB + 7) / 8 * 8 ))
if [ "$EST_GB" -lt 16 ];  then EST_GB=16;  fi
if [ "$EST_GB" -gt 300 ]; then EST_GB=300; fi
MEM="${EST_GB}G"

BASENAME="$(basename "$INPUT_FASTA")"
BASENAME="${BASENAME%.fasta}"; BASENAME="${BASENAME%.aln}"
BASENAME="${BASENAME%.fas}";   BASENAME="${BASENAME%.fa}"
BASENAME="${BASENAME%.fna}"

# ------------------------------------------------------------------------------
# 5. RESUMEN Y CONFIRMACIÓN
# ------------------------------------------------------------------------------
echo ""
linea
echo " RESUMEN DE LA CONFIGURACIÓN SELECCIONADA"
linea
echo " Secuencias:            $SEQ_TYPE_TXT ($NUM_SEQS secuencias, $ALN_LEN posiciones)"
echo " Modelo de sustitución:  $MODEL_TXT"
echo " Categorías de tasa:    $CAT_NUM categorías (-cat $CAT_NUM)"
echo " Optimización Gamma:    $GAMMA_TXT"
echo " Búsqueda topológica:   $SEARCH_TXT"
echo " Soporte de ramas:      $SUPPORT_TXT"
echo " Comand line flags:     FastTree $FASTTREE_EXTRA"
echo " Directorio salida:     $WORKDIR/${BASENAME}_fasttree_results/"
echo ""
echo " Reserva en Slurm:"
echo "   - Cómputo:  $CPUS CPUs (OpenMP FastTreeMP) | Memoria: $MEM"
echo "   - Tiempo:   máximo $TIME | Cola: $PARTITION"
if [ -n "$HUGE_NOTE" ]; then
    echo ""
    echo " $HUGE_NOTE"
fi
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
    --output="$WORKDIR/fasttree_%j.out" \
    --error="$WORKDIR/fasttree_%j.err" \
    --partition="$PARTITION" \
    --time="$TIME" \
    --cpus-per-task="$CPUS" \
    --mem="$MEM" \
    --job-name="FT_${NUM_SEQS}" \
    "$MASTER_SCRIPT" "$INPUT_FASTA" "$FASTTREE_EXTRA")
JOB_ID="${JOB_ID%%;*}"

echo ""
echo " Job enviado con éxito a Slurm. ID del trabajo: $JOB_ID"
echo ""
echo " Seguimiento y gestión:"
echo "   - Estado de la cola:            squeue -u \$USER"
echo "   - Ver salida en tiempo real:     tail -f $WORKDIR/fasttree_${JOB_ID}.err"
echo "   - Directorio de resultados:     $WORKDIR/${BASENAME}_fasttree_results/"