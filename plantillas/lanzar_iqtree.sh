#!/bin/bash
# ==============================================================================
# LANZADOR INTERACTIVO DE IQ-TREE 2 (VirResp-cluster)
# Uso: ./lanzar_iqtree.sh alineamiento.fasta
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
echo " ÁRBOL FILOGENÉTICO DE MÁXIMA VEROSIMILITUD CON IQ-TREE 2"
linea
echo " Fichero de entrada: Alineamiento FASTA (*.fasta, *.aln, *.fa, *.fas)"
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
MASTER_SCRIPT="$(realpath -m "$SCRIPT_DIR/../scripts/master_iqtree.sh")"

if [ ! -x "$MASTER_SCRIPT" ]; then
    error "No existe o no es ejecutable el script maestro:" "$MASTER_SCRIPT"
fi

# ------------------------------------------------------------------------------
# 2. INSPECCIÓN Y VALIDACIÓN DEL ALINEAMIENTO
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
    error "Se detectaron $NUM_SEQS secuencias." "Se requieren al menos 4 secuencias para construir un árbol filogenético."
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
# 3. SELECCIÓN DE PARÁMETROS DE EVOLUCIÓN Y MODELOS
# ------------------------------------------------------------------------------
RESIDUES=$(( TOTAL_CHARS - GAP_CHARS ))
if [ "$RESIDUES" -le 0 ]; then
    error "El alineamiento solo contiene huecos o caracteres ambiguos (-)."
fi
NT_PCT=$(( 100 * NUC_CHARS / RESIDUES ))

if [ "$NT_PCT" -ge 90 ]; then
    SEQ_TYPE="NT"
    DETECTED_TXT="Nucleótidos (ADN/ARN)"
else
    SEQ_TYPE="AA"
    DETECTED_TXT="Aminoácidos (Proteínas)"
fi

echo ""
echo "1. TIPO DE DATASET Y MODELO DE SUSTITUCIÓN EVOLUTIVA"
echo "   Autodetección basada en composición: $DETECTED_TXT"
echo ""

if [ "$SEQ_TYPE" = "NT" ]; then
    echo "  1) Búsqueda automática con ModelFinder (-m MFP) [Incluye FreeRate; Recomendado]"
    echo "  2) Búsqueda estándar con ModelFinder (-m MF) [Sin modelos FreeRate]"
    echo "  3) GTR+F+I+G4   - Modelo nucleotídico completo y complejo (ADN/ARN)"
    echo "  4) GTR+G        - Modelo estándar de Máxima Verosimilitud (ADN/ARN)"
    echo "  5) HKY85+G      - Modelo sencillo para secuencias poco divergentes"
    echo "  6) TN93+G       - Modelo de Tamura-Nei (frecuente en virus ARN)"
    echo "  7) Personalizado - Escribir modelo explícito de IQ-TREE"
    ask_choice "Selección" "1" 7
    MODEL_OPT="$CHOICE"

    case "$MODEL_OPT" in
        1) MODEL_FLAG="-m MFP";        MODEL_TXT="ModelFinder Plus (-m MFP)" ;;
        2) MODEL_FLAG="-m MF";         MODEL_TXT="ModelFinder Estándar (-m MF)" ;;
        3) MODEL_FLAG="-m GTR+F+I+G4"; MODEL_TXT="GTR+F+I+G4" ;;
        4) MODEL_FLAG="-m GTR+G";      MODEL_TXT="GTR+G" ;;
        5) MODEL_FLAG="-m HKY85+G";    MODEL_TXT="HKY85+G" ;;
        6) MODEL_FLAG="-m TN93+G";     MODEL_TXT="TN93+G" ;;
        7)
            read -r -p "   Introduce la sintaxis del modelo (ej. GTR+F+R4): " CUSTOM_MODEL
            MODEL_FLAG="-m ${CUSTOM_MODEL:-GTR+G}"
            MODEL_TXT="Personalizado ($MODEL_FLAG)"
            ;;
    esac
else
    echo "  1) Búsqueda automática con ModelFinder (-m MFP) [Recomendado para Proteínas]"
    echo "  2) Búsqueda estándar con ModelFinder (-m MF)"
    echo "  3) LG+F+G4      - Modelo estándar para aminoácidos virales"
    echo "  4) WAG+F+G4     - Modelo Whelan & Goldman"
    echo "  5) JTT+F+G4     - Modelo Jones-Taylor-Thornton"
    echo "  6) Personalizado - Escribir modelo explícito de IQ-TREE"
    ask_choice "Selección" "1" 6
    MODEL_OPT="$CHOICE"

    case "$MODEL_OPT" in
        1) MODEL_FLAG="-m MFP";       MODEL_TXT="ModelFinder Plus (-m MFP)" ;;
        2) MODEL_FLAG="-m MF";        MODEL_TXT="ModelFinder Estándar (-m MF)" ;;
        3) MODEL_FLAG="-m LG+F+G4";   MODEL_TXT="LG+F+G4" ;;
        4) MODEL_FLAG="-m WAG+F+G4";  MODEL_TXT="WAG+F+G4" ;;
        5) MODEL_FLAG="-m JTT+F+G4";  MODEL_TXT="JTT+F+G4" ;;
        6)
            read -r -p "   Introduce la sintaxis del modelo (ej. LG+C20+F+G): " CUSTOM_MODEL
            MODEL_FLAG="-m ${CUSTOM_MODEL:-LG+F+G4}"
            MODEL_TXT="Personalizado ($MODEL_FLAG)"
            ;;
    esac
fi

echo ""
echo "2. SOPORTE DE RAMAS Y MÉTODOS DE EVALUACIÓN (BOOTSTRAP)"
echo "  1) Ultrafast Bootstrap (UFBoot) - 1000 réplicas (-bb 1000) [Rápido y Estándar]"
echo "  2) Doble test: SH-aLRT + UFBoot (-alrt 1000 -bb 1000) [Máxima Rigurosidad]"
echo "  3) Test SH-aLRT únicamente (-alrt 1000) [Rápido en datasets grandes]"
echo "  4) Bootstrap Estándar No Paramétrico (-b 100) [Muy lento; solo conjuntos pequeños]"
echo "  5) Sin soporte de ramas (generar únicamente la topología ML)"
ask_choice "Selección" "1" 5
BOOT_OPT="$CHOICE"

case "$BOOT_OPT" in
    1) BOOT_FLAGS="-bb 1000 -bnni";            BOOT_TXT="UFBoot 1000 réplicas + Optimización BNNI (-bb 1000 -bnni)" ;;
    2) BOOT_FLAGS="-alrt 1000 -bb 1000 -bnni"; BOOT_TXT="SH-aLRT (1000) + UFBoot (1000) + BNNI" ;;
    3) BOOT_FLAGS="-alrt 1000";                 BOOT_TXT="SH-aLRT (1000 réplicas)" ;;
    4)
        if [ "$NUM_SEQS" -gt 500 ]; then
            echo ""
            echo " [AVISO] El Bootstrap estándar (-b 100) es excesivamente lento para $NUM_SEQS secuencias."
            if ! ask_yn "         ¿Deseas cambiar a Ultrafast Bootstrap (-bb 1000)?" "s"; then
                BOOT_FLAGS="-b 100"
                BOOT_TXT="Bootstrap Estándar (100 réplicas)"
            else
                BOOT_FLAGS="-bb 1000 -bnni"
                BOOT_TXT="UFBoot 1000 réplicas + BNNI"
            fi
        else
            BOOT_FLAGS="-b 100"
            BOOT_TXT="Bootstrap Estándar (100 réplicas)"
        fi
        ;;
    5) BOOT_FLAGS=""; BOOT_TXT="Sin soporte de ramas" ;;
esac

echo ""
echo "3. FIJACIÓN DE SEMILLA ALEATORIA PARA REPRODUCIBILIDAD (--seed)"
if ask_yn "   ¿Utilizar una semilla aleatoria fija (seed 1253)?" "s"; then
    SEED_FLAG="-seed 1253"
    SEED_TXT="Semilla fija (1253)"
else
    SEED_FLAG=""
    SEED_TXT="Semilla aleatoria del sistema"
fi

# ------------------------------------------------------------------------------
# 4. CONSTRUCCIÓN DE COMANDOS Y ESTIMACIÓN DE RECURSOS SLURM
# ------------------------------------------------------------------------------
FLAGS=()
FLAGS+=($MODEL_FLAG)

if [ -n "$BOOT_FLAGS" ]; then
    FLAGS+=($BOOT_FLAGS)
fi

if [ -n "$SEED_FLAG" ]; then
    FLAGS+=($SEED_FLAG)
fi

# Factor de esfuerzo para Slurm
WORK=$(( NUM_SEQS * ALN_LEN ))
IS_MF=0
if [[ "$MODEL_FLAG" == *"-m MF"* ]]; then IS_MF=1; fi

if [ "$NUM_SEQS" -gt 5000 ] || [ "$WORK" -gt 2000000000 ]; then
    PARTITION="long_idx";   TIME="5-00:00:00"; CPUS=32; MEM_GB=300
elif [ "$NUM_SEQS" -gt 1000 ] || [ "$IS_MF" -eq 1 ] || [ "$WORK" -gt 100000000 ]; then
    PARTITION="middle_idx"; TIME="48:00:00";   CPUS=16; MEM_GB=128
else
    PARTITION="short_idx";  TIME="12:00:00";   CPUS=16; MEM_GB=64
fi

# Añadir límite de memoria explícito a IQ-TREE 2
FLAGS+=("-mem" "${MEM_GB}G")
IQTREE_EXTRA="${FLAGS[*]}"
MEM="${MEM_GB}G"

BASENAME="$(basename "$INPUT_FASTA")"
BASENAME="${BASENAME%.fasta}"; BASENAME="${BASENAME%.aln}"
BASENAME="${BASENAME%.fas}";   BASENAME="${BASENAME%.fa}"
BASENAME="${BASENAME%.fna}";   BASENAME="${BASENAME%.aligned}"

# ------------------------------------------------------------------------------
# 5. RESUMEN Y CONFIRMACIÓN
# ------------------------------------------------------------------------------
echo ""
linea
echo " RESUMEN DE LA CONFIGURACIÓN DE IQ-TREE 2"
linea
echo " Archivo de entrada:     $INPUT_FASTA ($NUM_SEQS secuencias, $ALN_LEN pos)"
echo " Tipo de datos:          $DETECTED_TXT"
echo " Modelo evolutivo:       $MODEL_TXT"
echo " Soporte de ramas:       $BOOT_TXT"
echo " Reproducibilidad:       $SEED_TXT"
echo " Comand line flags:      iqtree2 -s input.fasta -nt AUTO -ntmax \$SLURM_CPUS_PER_TASK $IQTREE_EXTRA"
echo " Directorio salida:      $WORKDIR/${BASENAME}_iqtree_results/"
echo ""
echo " Reserva en Slurm:"
echo "   - Cómputo:  $CPUS CPUs (multithread OpenMP) | Memoria RAM: $MEM"
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
    --output="$WORKDIR/iqtree_%j.out" \
    --error="$WORKDIR/iqtree_%j.err" \
    --partition="$PARTITION" \
    --time="$TIME" \
    --cpus-per-task="$CPUS" \
    --mem="$MEM" \
    --job-name="IQT_${NUM_SEQS}" \
    "$MASTER_SCRIPT" "$INPUT_FASTA" "$IQTREE_EXTRA")
JOB_ID="${JOB_ID%%;*}"

echo ""
echo " Job enviado con éxito a Slurm. ID del trabajo: $JOB_ID"
echo ""
echo " Seguimiento y gestión:"
echo "   - Estado de la cola:            squeue -u \$USER"
echo "   - Ver salida en tiempo real:     tail -f $WORKDIR/iqtree_${JOB_ID}.err"
echo "   - Directorio de resultados:     $WORKDIR/${BASENAME}_iqtree_results/"