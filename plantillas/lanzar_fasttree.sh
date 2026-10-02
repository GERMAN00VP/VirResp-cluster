#!/bin/bash
# ==============================================================================
# LANZADOR INTERACTIVO DE FASTTREE  (VirResp-cluster)
# Uso:  ./lanzar_fasttree.sh alineamiento.fasta
#
# Qué hace: revisa tu alineamiento, te hace 3 preguntas sencillas, calcula
# solo los recursos que hacen falta y envía el trabajo al clúster.
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

# Pregunta con opciones numéricas. Resultado en la variable CHOICE.
#   $1 = texto, $2 = opción por defecto, $3 = número máximo de opciones
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

# Pregunta sí/no. Devuelve 0 si la respuesta es sí.
#   $1 = texto, $2 = respuesta por defecto (s o n)
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
# Mensaje inicial
# ------------------------------------------------------------------------------
linea
echo " ÁRBOL FILOGENÉTICO RÁPIDO CON FASTTREE"
linea
echo " Necesitas un archivo FASTA ya ALINEADO (por ejemplo, con MAFFT)."
echo " Consejo: ejecuta este script desde la carpeta donde está tu archivo."
echo " Los resultados aparecerán en una carpeta llamada *_fasttree_results/"
linea
echo ""

# ------------------------------------------------------------------------------
# 1. COMPROBACIONES PREVIAS
# ------------------------------------------------------------------------------
if [ "$#" -ne 1 ]; then
    error "Tienes que indicar el archivo alineado." "Uso: $0 <alineamiento.fasta>"
fi

if [ ! -f "$1" ]; then
    error "El archivo '$1' no existe." "Comprueba el nombre con:  ls -l"
fi

if [ ! -s "$1" ]; then
    error "El archivo '$1' está vacío."
fi

case "$1" in
    *.gz) error "El archivo está comprimido (.gz)." "Descomprímelo primero con:  gunzip $1" ;;
esac

if ! command -v sbatch >/dev/null 2>&1; then
    error "No encuentro el gestor de trabajos del clúster (sbatch)." \
          "Este script solo funciona conectado al clúster."
fi

INPUT_FASTA="$(realpath "$1")"
WORKDIR="$(dirname "$INPUT_FASTA")"

if [ ! -w "$WORKDIR" ]; then
    error "No tienes permiso para escribir en la carpeta:" "$WORKDIR" \
          "Copia el archivo a tu carpeta de usuario y vuelve a probar."
fi

# El script maestro vive junto a este lanzador (../scripts/), esté donde esté el repo
SCRIPT_DIR="$(cd "$(dirname "$(realpath "${BASH_SOURCE[0]}")")" && pwd)"
MASTER_SCRIPT="$(realpath -m "$SCRIPT_DIR/../scripts/master_fasttree.sh")"

if [ ! -x "$MASTER_SCRIPT" ]; then
    error "Hay un problema con la instalación del repositorio (no es ejecutable):" \
          "$MASTER_SCRIPT" "Avisa al administrador (Germán)."
fi

# ------------------------------------------------------------------------------
# 2. ANÁLISIS DEL ALINEAMIENTO (no modifica tu archivo)
# ------------------------------------------------------------------------------
# Una sola lectura: nº de secuencias, longitud mínima y máxima, y qué tipo de
# letras contiene (para adivinar si son nucleótidos o proteínas).
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
    error "El archivo no parece un FASTA: no encuentro ninguna secuencia (líneas que empiezan por '>')."
fi

if [ "$NUM_SEQS" -lt 4 ]; then
    error "Hay solo $NUM_SEQS secuencias." "Para construir un árbol se necesitan al menos 4."
fi

ALN_LEN="$MAX_LEN"

echo " ANALIZANDO TU ALINEAMIENTO"
echo " -> Archivo:      $INPUT_FASTA"
echo " -> Secuencias:   $NUM_SEQS"
echo " -> Longitud:     $ALN_LEN posiciones"

# Todas las secuencias deben medir lo mismo si están alineadas
if [ "$MIN_LEN" -ne "$MAX_LEN" ]; then
    error "Tus secuencias NO miden todas lo mismo (de $MIN_LEN a $MAX_LEN posiciones)." \
          "Eso significa que el archivo NO está alineado." \
          "Primero alinéalo con:  lanzar_mafft.sh  y usa el archivo alineado que genera."
fi

# Nombres repetidos
DUPS=$(tr -d '\r' < "$INPUT_FASTA" | grep '^>' | sed 's/^>//; s/[[:space:]].*$//' | sort | uniq -d | head -n 5 || true)
if [ -n "$DUPS" ]; then
    error "Hay secuencias con el MISMO nombre. Cada secuencia necesita un nombre único." \
          "Algunos ejemplos repetidos: $(echo "$DUPS" | tr '\n' ' ')"
fi

# Caracteres que dan problemas en los árboles
BAD_IDS=$(tr -d '\r' < "$INPUT_FASTA" | grep -c -E '^>[^[:space:]]*[][(),:;]' || true)
if [ "$BAD_IDS" -gt 0 ]; then
    echo ""
    echo " [AVISO] $BAD_IDS nombres contienen alguno de estos símbolos:  ( ) , : ; [ ]"
    echo "         Pueden estropear el archivo del árbol. Te recomiendo cambiarlos"
    echo "         (por ejemplo por guiones bajos) y volver a empezar."
    if ! ask_yn "         ¿Quieres continuar de todos modos?" "n"; then
        echo "Operación cancelada."
        exit 0
    fi
fi
linea

# ------------------------------------------------------------------------------
# 3. PREGUNTA 1: TIPO DE SECUENCIA
# ------------------------------------------------------------------------------
RESIDUES=$(( TOTAL_CHARS - GAP_CHARS ))
if [ "$RESIDUES" -le 0 ]; then
    error "El alineamiento solo contiene huecos (-). No hay nada que analizar."
fi
NT_PCT=$(( 100 * NUC_CHARS / RESIDUES ))

if [ "$NT_PCT" -ge 90 ]; then
    DETECTED=1
    DETECTED_TXT="nucleótidos (ADN/ARN)"
else
    DETECTED=2
    DETECTED_TXT="aminoácidos (proteínas)"
fi

echo ""
echo "1. ¿Qué tipo de secuencias son?"
echo "   Por su contenido, parecen ser: $DETECTED_TXT"
echo "  1) Nucleótidos (ADN / ARN): genomas o genes de virus"
echo "  2) Aminoácidos (proteínas)"
ask_choice "Opción" "$DETECTED" 2
SEQ_OPT="$CHOICE"

if [ "$SEQ_OPT" -ne "$DETECTED" ]; then
    echo ""
    echo " [AVISO] Tu archivo parece contener $DETECTED_TXT, pero has elegido otra cosa."
    if ! ask_yn "         ¿Seguro que quieres continuar con tu elección?" "n"; then
        echo "Operación cancelada."
        exit 0
    fi
fi

if [ "$SEQ_OPT" -eq 1 ]; then
    SEQ_TYPE="NT"
    SEQ_TYPE_TXT="Nucleótidos (ADN/ARN)"
    MODEL_TXT="GTR (el modelo estándar para virus de ADN/ARN; se aplica solo)"
else
    SEQ_TYPE="AA"
    SEQ_TYPE_TXT="Aminoácidos (proteínas)"
    MODEL_TXT="LG (el modelo estándar para proteínas; se aplica solo)"
fi

# ------------------------------------------------------------------------------
# 4. PREGUNTA 2: CUIDADO EN LA BÚSQUEDA DEL ÁRBOL
# ------------------------------------------------------------------------------
echo ""
echo "2. ¿Qué tipo de búsqueda quieres?"
MODE_OPT=1
if [ "$NUM_SEQS" -le 20000 ]; then
    echo "  1) Estándar (recomendado): rápido y fiable para casi todos los casos."
    echo "  2) Cuidadoso: el programa dedica más tiempo a buscar el mejor árbol posible"
    echo "     (tarda varias veces más). Útil para el árbol definitivo de un artículo."
    ask_choice "Opción" "1" 2
    MODE_OPT="$CHOICE"
else
    echo "   Con tantas secuencias ($NUM_SEQS) se usa automáticamente la búsqueda estándar."
fi

# ------------------------------------------------------------------------------
# 5. PREGUNTA 3: FIABILIDAD DE LAS RAMAS
# ------------------------------------------------------------------------------
echo ""
echo "3. ¿Cómo quieres medir la fiabilidad de cada rama del árbol?"
echo "   (Son los números que aparecen sobre las ramas, de 0 a 1.)"
if [ "$NUM_SEQS" -le 5000 ]; then
    echo "  1) Rápida (recomendado): valores 'SH-like', se calculan sin coste extra."
    echo "     No equivalen a un bootstrap, así que no los llames 'bootstrap' en un artículo."
    echo "  2) Bootstrap (100 repeticiones): el método clásico que piden muchas revistas."
    echo "     Tarda bastante más."
    echo "  3) Sin valores de fiabilidad: lo más rápido."
    ask_choice "Opción" "1" 3
    SUPPORT_OPT="$CHOICE"
else
    echo "  1) Rápida (recomendado): valores 'SH-like', se calculan sin coste extra."
    echo "     No equivalen a un bootstrap, así que no los llames 'bootstrap' en un artículo."
    echo "  2) Sin valores de fiabilidad: lo más rápido."
    echo "   (El bootstrap no está disponible con más de 5.000 secuencias.)"
    ask_choice "Opción" "1" 2
    if [ "$CHOICE" -eq 2 ]; then SUPPORT_OPT=3; else SUPPORT_OPT=1; fi
fi

case "$SUPPORT_OPT" in
    1) SUPPORT_TXT="Rápida (SH-like)" ;;
    2) SUPPORT_TXT="Bootstrap (100 repeticiones)" ;;
    3) SUPPORT_TXT="Sin valores de fiabilidad" ;;
esac
if [ "$MODE_OPT" -eq 2 ]; then MODE_TXT="Cuidadoso"; else MODE_TXT="Estándar"; fi

# ------------------------------------------------------------------------------
# 6. CONFIGURACIÓN TÉCNICA (automática; el usuario no tiene que decidir nada)
# ------------------------------------------------------------------------------
FLAGS=()
if [ "$SEQ_TYPE" = "NT" ]; then
    FLAGS+=(-nt -gtr)
else
    FLAGS+=(-lg)
fi
FLAGS+=(-gamma -seed 1253)

HUGE_NOTE=""
if [ "$NUM_SEQS" -gt 50000 ]; then
    FLAGS+=(-fastest)
    HUGE_NOTE="Conjunto muy grande (>50.000 secuencias): se activa un modo de ahorro de memoria y tiempo."
fi
if [ "$MODE_OPT" -eq 2 ]; then
    FLAGS+=(-spr 4 -mlacc 2 -slownni)
fi
case "$SUPPORT_OPT" in
    2) FLAGS+=(-boot 100) ;;
    3) FLAGS+=(-nosupport) ;;
esac
FASTTREE_EXTRA="${FLAGS[*]}"

# Estimación de recursos según secuencias x longitud (no solo nº de secuencias)
WORK=$(( NUM_SEQS * ALN_LEN ))
EFF_WORK="$WORK"
if [ "$MODE_OPT" -eq 2 ];    then EFF_WORK=$(( EFF_WORK * 3 )); fi
if [ "$SUPPORT_OPT" -eq 2 ]; then EFF_WORK=$(( EFF_WORK * 4 )); fi

if [ "$NUM_SEQS" -gt 50000 ] || [ "$EFF_WORK" -gt 1500000000 ]; then
    PARTITION="long_idx";   TIME="5-00:00:00"; CPUS=16
elif [ "$NUM_SEQS" -gt 5000 ] || [ "$EFF_WORK" -gt 100000000 ]; then
    PARTITION="middle_idx"; TIME="48:00:00";   CPUS=8
else
    PARTITION="short_idx";  TIME="12:00:00";   CPUS=4
fi

# Memoria orientativa: 2 nodos por secuencia, 16 bytes/posición (NT) o 80 (AA),
# x1,5 de margen + 5 GB base. Se redondea a múltiplos de 8 GB (mín. 16, máx. 300).
if [ "$SEQ_TYPE" = "NT" ]; then BYTES_POS=16; else BYTES_POS=80; fi
EST_BYTES=$(( 3 * NUM_SEQS * ALN_LEN * BYTES_POS ))
EST_GB=$(( EST_BYTES / 1000000000 + 5 ))
EST_GB=$(( (EST_GB + 7) / 8 * 8 ))
if [ "$EST_GB" -lt 16 ];  then EST_GB=16;  fi
if [ "$EST_GB" -gt 300 ]; then EST_GB=300; fi
MEM="${EST_GB}G"

# Nombre base (igual que en el script maestro) para anunciar la carpeta de resultados
BASENAME="$(basename "$INPUT_FASTA")"
BASENAME="${BASENAME%.fasta}"; BASENAME="${BASENAME%.fas}"
BASENAME="${BASENAME%.fa}";    BASENAME="${BASENAME%.fna}"

# ------------------------------------------------------------------------------
# 7. RESUMEN Y CONFIRMACIÓN
# ------------------------------------------------------------------------------
echo ""
linea
echo " RESUMEN"
linea
echo " Secuencias:           $SEQ_TYPE_TXT"
echo " Modelo evolutivo:     $MODEL_TXT"
echo " Tipo de búsqueda:     $MODE_TXT"
echo " Fiabilidad de ramas:  $SUPPORT_TXT"
echo " Carpeta de salida:    $WORKDIR/${BASENAME}_fasttree_results/"
echo ""
echo " Se reservará en el clúster:"
echo "   - $CPUS núcleos de cálculo y $MEM de memoria"
echo "   - hasta $TIME de tiempo (cola '$PARTITION')"
echo "   (Es una estimación a partir del tamaño de tu alineamiento.)"
if [ -n "$HUGE_NOTE" ]; then
    echo ""
    echo " $HUGE_NOTE"
fi
linea

if ! ask_yn "¿Enviar el trabajo al clúster?" "s"; then
    echo "Operación cancelada. No se ha enviado nada."
    exit 0
fi

# ------------------------------------------------------------------------------
# 8. ENVÍO AL CLÚSTER
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
echo " Trabajo enviado correctamente. Número de trabajo: $JOB_ID"
echo ""
echo " Qué hacer ahora:"
echo "   - Ver si ya está en marcha o terminado:   squeue -u \$USER"
echo "     (si no aparece, ha terminado: mira la carpeta de resultados)"
echo "   - Seguir lo que va haciendo:              tail -f $WORKDIR/fasttree_${JOB_ID}.err"
echo "     (ese archivo .err muestra mensajes de progreso normales; no significa que haya fallos)"
echo "   - Resultados al terminar:                 $WORKDIR/${BASENAME}_fasttree_results/"
