#!/bin/bash
# ==============================================================================
# LANZADOR INTERACTIVO DE IQ-TREE 2 (VirResp-cluster)
# Uso: ./lanzar_iqtree.sh alineamiento.fasta
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

MASTER_SCRIPT="$(realpath -m "$SCRIPT_DIR/../scripts/master_iqtree.sh")"

if [ ! -x "$MASTER_SCRIPT" ]; then
    error "No existe o no es ejecutable el script maestro:" "$MASTER_SCRIPT"
fi

# ------------------------------------------------------------------------------
# 2. INSPECCIÓN Y VALIDACIÓN DEL ALINEAMIENTO
# ------------------------------------------------------------------------------
fasta_basico "$INPUT_FASTA"

if [ "$NUM_SEQS" -eq 0 ]; then
    error "No se detectaron cabeceras FASTA (líneas '>') en '$INPUT_FASTA'."
fi

if [ "$NUM_SEQS" -lt 4 ]; then
    error "Se detectaron $NUM_SEQS secuencias." "Se requieren al menos 4 secuencias para construir un árbol filogenético."
fi

fasta_longitudes_muestra "$INPUT_FASTA" 50

echo " INSPECCIÓN DEL ALINEAMIENTO"
echo " -> Archivo:      $INPUT_FASTA"
echo " -> Secuencias:   $NUM_SEQS"
echo " -> Longitud:     $MAX_LEN pb/aa (comprobada en $N_MUESTRA primeras + última)"

fasta_exigir_alineado

DUPS=$(fasta_ids_duplicados "$INPUT_FASTA")
if [ -n "$DUPS" ]; then
    error "Identificadores de secuencia duplicados detectados:" "$(echo "$DUPS" | tr '\n' ' ')"
fi

BAD_IDS=$(fasta_ids_newick "$INPUT_FASTA")
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
fasta_tipo_muestra "$INPUT_FASTA"

if [ "$NT_PCT" -ge 90 ]; then
    SEQ_TYPE="NT"
    DETECTED_TXT="Nucleótidos (ADN/ARN)"
else
    SEQ_TYPE="AA"
    DETECTED_TXT="Aminoácidos (Proteínas)"
fi

echo ""
echo "1. TIPO DE DATASET Y MODELO DE SUSTITUCIÓN EVOLUTIVA"
echo "   Autodetección basada en composición (muestra de ~2 MB): $DETECTED_TXT"
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

# Aviso de ModelFinder en datasets grandes
IS_MF=0
if [[ "$MODEL_FLAG" == *"-m MF"* ]]; then IS_MF=1; fi
if [ "$IS_MF" -eq 1 ] && [ "$NUM_SEQS" -gt 3000 ]; then
    if [ "$SEQ_TYPE" = "NT" ]; then ALT_MODEL="GTR+F+R4"; else ALT_MODEL="LG+F+G4"; fi
    echo ""
    echo " [AVISO] ModelFinder evalúa decenas de modelos y cada uno reoptimiza el árbol completo."
    echo "         Con ~18 800 secuencias de SARS-CoV-2 tardó entre 1 y 5 h por modelo"
    echo "         (7 modelos en 16 h). Alternativas: un modelo fijo (GTR+F+R4 fue el mejor"
    echo "         por BIC en ese caso real; también GTR+F+I+G4) o ejecutar ModelFinder antes"
    echo "         sobre un subconjunto de ~500 secuencias y usar aquí el modelo elegido."
    if ask_yn "         ¿Cambiar al modelo fijo $ALT_MODEL?" "s"; then
        MODEL_FLAG="-m $ALT_MODEL"
        MODEL_TXT="$ALT_MODEL (modelo fijo en lugar de ModelFinder)"
        IS_MF=0
    fi
fi

echo ""
echo "2. TIPO DE BÚSQUEDA DEL ÁRBOL"
if [ "$NUM_SEQS" -gt 10000 ]; then
    echo " [AVISO] Dataset de escala pandémica ($NUM_SEQS secuencias). IQ-TREE puede tardar"
    echo "         días incluso en modo rápido. Para este tamaño se recomiendan herramientas"
    echo "         específicas como MAPLE, UShER o VeryFastTree."
fi
echo "  1) Estándar (100 árboles iniciales + búsqueda estocástica)"
echo "  2) Rápida (-fast) [Recomendado si > 3000 secuencias]"
if [ "$NUM_SEQS" -gt 3000 ]; then DEF_SEARCH=2; else DEF_SEARCH=1; fi
ask_choice "Selección" "$DEF_SEARCH" 2
FAST=0
if [ "$CHOICE" -eq 2 ]; then FAST=1; fi

if [ "$FAST" -eq 0 ] && [ "$NUM_SEQS" -gt 3000 ]; then
    echo ""
    echo " [AVISO] La búsqueda estándar genera 100 árboles de parsimonia (en un solo núcleo;"
    echo "         ~20 min cada uno con 18 800 secuencias) y luego optimiza los candidatos."
    echo "         Puede superar los 5 días del límite de cola."
    if ask_yn "         ¿Cambiar a búsqueda rápida (-fast)?" "s"; then
        FAST=1
    fi
fi

echo ""
echo "3. SOPORTE DE RAMAS Y MÉTODOS DE EVALUACIÓN (BOOTSTRAP)"
if [ "$FAST" -eq 1 ]; then
    echo "   Con -fast, IQ-TREE no admite UFBoot (-bb); se recomienda SH-aLRT."
    echo "  1) Ultrafast Bootstrap (UFBoot) - 1000 réplicas (-bb 1000) [No compatible con -fast]"
    echo "  2) Doble test: SH-aLRT + UFBoot (-alrt 1000 -bb 1000) [No compatible con -fast]"
    echo "  3) Test SH-aLRT únicamente (-alrt 1000) [Recomendado con -fast]"
    DEF_BOOT=3
else
    echo "  1) Ultrafast Bootstrap (UFBoot) - 1000 réplicas (-bb 1000) [Rápido y Estándar]"
    echo "  2) Doble test: SH-aLRT + UFBoot (-alrt 1000 -bb 1000) [Máxima Rigurosidad]"
    echo "  3) Test SH-aLRT únicamente (-alrt 1000) [Rápido en datasets grandes]"
    DEF_BOOT=1
fi
echo "  4) Bootstrap Estándar No Paramétrico (-b 100) [Muy lento; solo conjuntos pequeños]"
echo "  5) Sin soporte de ramas (generar únicamente la topología ML)"
ask_choice "Selección" "$DEF_BOOT" 5
BOOT_OPT="$CHOICE"

# UFBoot y -fast son incompatibles en IQ-TREE ("Ultrafast bootstrap does not work with -fast")
if [ "$FAST" -eq 1 ] && { [ "$BOOT_OPT" -eq 1 ] || [ "$BOOT_OPT" -eq 2 ]; }; then
    echo ""
    echo " [AVISO] UFBoot no funciona junto con la búsqueda rápida (-fast)."
    if ask_yn "         ¿Usar SH-aLRT (-alrt 1000) y mantener -fast? (no = búsqueda estándar)" "s"; then
        BOOT_OPT=3
    else
        FAST=0
    fi
fi

case "$BOOT_OPT" in
    1) BOOT_FLAGS="-bb 1000 -bnni";            BOOT_TXT="UFBoot 1000 réplicas + Optimización BNNI (-bb 1000 -bnni)" ;;
    2) BOOT_FLAGS="-alrt 1000 -bb 1000 -bnni"; BOOT_TXT="SH-aLRT (1000) + UFBoot (1000) + BNNI" ;;
    3) BOOT_FLAGS="-alrt 1000";                 BOOT_TXT="SH-aLRT (1000 réplicas)" ;;
    4)
        if [ "$NUM_SEQS" -gt 500 ]; then
            echo ""
            echo " [AVISO] El Bootstrap estándar (-b 100) es excesivamente lento para $NUM_SEQS secuencias."
            if [ "$FAST" -eq 1 ]; then
                ALT_BOOT="-alrt 1000"; ALT_BOOT_TXT="SH-aLRT (1000 réplicas)"
            else
                ALT_BOOT="-bb 1000 -bnni"; ALT_BOOT_TXT="UFBoot 1000 réplicas + BNNI"
            fi
            if ! ask_yn "         ¿Deseas cambiar a $ALT_BOOT_TXT ($ALT_BOOT)?" "s"; then
                BOOT_FLAGS="-b 100"
                BOOT_TXT="Bootstrap Estándar (100 réplicas)"
            else
                BOOT_FLAGS="$ALT_BOOT"
                BOOT_TXT="$ALT_BOOT_TXT"
            fi
        else
            BOOT_FLAGS="-b 100"
            BOOT_TXT="Bootstrap Estándar (100 réplicas)"
        fi
        ;;
    5) BOOT_FLAGS=""; BOOT_TXT="Sin soporte de ramas" ;;
esac

if [[ " $BOOT_FLAGS " == *" -bb "* ]] && [ "$NUM_SEQS" -gt 3000 ]; then
    echo ""
    echo " [AVISO] Con muchos taxones ($NUM_SEQS) la búsqueda con UFBoot puede no converger"
    echo "         dentro del límite de la cola."
    if ask_yn "         ¿Cambiar a SH-aLRT únicamente (-alrt 1000)?" "s"; then
        BOOT_FLAGS="-alrt 1000"
        BOOT_TXT="SH-aLRT (1000 réplicas)"
    fi
fi

if [[ " $BOOT_FLAGS " == *" -bnni "* ]] && [ "$NUM_SEQS" -gt 1000 ]; then
    echo ""
    echo " [AVISO] -bnni añade una optimización NNI por cada réplica de UFBoot; es muy"
    echo "         costoso con $NUM_SEQS secuencias."
    if ask_yn "         ¿Quitar -bnni?" "s"; then
        BOOT_FLAGS="${BOOT_FLAGS% -bnni}"
        BOOT_TXT="$BOOT_TXT (sin BNNI)"
    fi
fi

if [ "$FAST" -eq 1 ]; then
    SEARCH_FLAG="-fast"
    SEARCH_TXT="Rápida (-fast)"
else
    SEARCH_FLAG=""
    SEARCH_TXT="Estándar (100 árboles iniciales + búsqueda estocástica)"
fi

if [ "$NUM_SEQS" -gt 200 ]; then DEF_THR=2; else DEF_THR=1; fi
echo ""
echo "4. HILOS DE CÓMPUTO"
echo "   IQ-TREE reparte el trabajo por patrones de sitio; en alineamientos con poca"
echo "   diversidad (p. ej. SARS-CoV-2) más de 4-8 hilos no acelera y solo bloquea CPUs."
echo "   Algunas fases (los árboles de parsimonia iniciales) usan un solo núcleo, y"
echo "   reservar más CPUs solo alarga la espera en la cola."
if [ "$DEF_THR" -eq 2 ]; then
    echo "  1) 4 hilos"
    echo "  2) 8 hilos   [Recomendado]"
else
    echo "  1) 4 hilos   [Recomendado]"
    echo "  2) 8 hilos"
fi
echo "  3) 16 hilos"
echo "  4) AUTO (IQ-TREE mide la eficiencia, máximo 8) [añade ~20-30 min de test en datasets grandes]"
echo "  5) Personalizado (1-32)"
ask_choice "Selección" "$DEF_THR" 5

THREAD_FLAGS=""
case "$CHOICE" in
    1) CPUS=4 ;;
    2) CPUS=8 ;;
    3) CPUS=16 ;;
    4) CPUS=8; THREAD_FLAGS="-T AUTO --threads-max 8" ;;
    5) pedir_numero_hilos; CPUS="$HILOS" ;;
esac
if [ -n "$THREAD_FLAGS" ]; then
    THREAD_TXT="AUTO (máximo 8 hilos; -T AUTO --threads-max 8)"
else
    THREAD_TXT="$CPUS hilos (-T $CPUS)"
fi

echo ""
echo "5. FIJACIÓN DE SEMILLA ALEATORIA PARA REPRODUCIBILIDAD (--seed)"
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

if [ -n "$SEARCH_FLAG" ]; then
    FLAGS+=("$SEARCH_FLAG")
fi

if [ -n "$BOOT_FLAGS" ]; then
    FLAGS+=($BOOT_FLAGS)
fi

if [ -n "$SEED_FLAG" ]; then
    FLAGS+=($SEED_FLAG)
fi

# Con AUTO el lanzador pasa -T y el maestro no lo añade
if [ -n "$THREAD_FLAGS" ]; then
    FLAGS+=($THREAD_FLAGS)
fi

# Estimación orientativa del tiempo (0 minutos, 1 horas, 2 días, 3 probablemente no termina)
HEAVY=0
if [ "$IS_MF" -eq 1 ] || [ "$FAST" -eq 0 ] || [[ " $BOOT_FLAGS " == *" -bb "* || " $BOOT_FLAGS " == *" -b "* ]]; then
    HEAVY=1
fi
if [ "$NUM_SEQS" -lt 500 ]; then
    EST_CODE=0
elif [ "$NUM_SEQS" -le 3000 ]; then
    EST_CODE=1
elif [ "$NUM_SEQS" -le 10000 ]; then
    if [ "$HEAVY" -eq 1 ]; then EST_CODE=2; else EST_CODE=1; fi
else
    if [ "$IS_MF" -eq 1 ] || [ "$FAST" -eq 0 ]; then EST_CODE=3; else EST_CODE=2; fi
fi
# El bootstrap no paramétrico multiplica el tiempo ~100 veces
if [[ " $BOOT_FLAGS " == *" -b "* ]] && [ "$EST_CODE" -lt 2 ]; then
    EST_CODE=$(( EST_CODE + 1 ))
fi
case "$EST_CODE" in
    0) EST_TXT="Minutos" ;;
    1) EST_TXT="Horas" ;;
    2) EST_TXT="1-3 días" ;;
    3) EST_TXT="Probablemente no termina en 5 días" ;;
esac

# Cola y tiempo por tramos; si se estiman más de 2 días, cola larga
WORK=$(( NUM_SEQS * ALN_LEN ))
if [ "$NUM_SEQS" -gt 5000 ] || [ "$WORK" -gt 2000000000 ] || [ "$EST_CODE" -ge 2 ]; then
    PARTITION="long_idx";   TIME="5-00:00:00"
elif [ "$NUM_SEQS" -gt 1000 ] || [ "$IS_MF" -eq 1 ] || [ "$WORK" -gt 100000000 ]; then
    PARTITION="middle_idx"; TIME="48:00:00"
else
    PARTITION="short_idx";  TIME="12:00:00"
fi

# Memoria por tamaño del alineamiento (independiente de los hilos)
if [ "$WORK" -le 100000000 ]; then
    MEM_GB=64
elif [ "$WORK" -le 1000000000 ]; then
    MEM_GB=128
else
    MEM_GB=256
fi

# Límite de memoria de IQ-TREE ~8G por debajo de lo reservado en Slurm
FLAGS+=("-mem" "$(( MEM_GB - 8 ))G")
IQTREE_EXTRA="${FLAGS[*]}"
MEM="${MEM_GB}G"

BASENAME="$(basename "$INPUT_FASTA")"
BASENAME="${BASENAME%.fasta}"; BASENAME="${BASENAME%.aln}"
BASENAME="${BASENAME%.fas}";   BASENAME="${BASENAME%.fa}"
BASENAME="${BASENAME%.fna}";   BASENAME="${BASENAME%.aligned}"

pedir_outdir "$WORKDIR/${BASENAME}_iqtree_results" "$WORKDIR"

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
echo " Búsqueda del árbol:     $SEARCH_TXT"
echo " Soporte de ramas:       $BOOT_TXT"
echo " Reproducibilidad:       $SEED_TXT"
echo " Hilos:                  $THREAD_TXT"
if [ -n "$THREAD_FLAGS" ]; then
    echo " Comand line flags:      iqtree2 -s input.fasta $IQTREE_EXTRA"
else
    echo " Comand line flags:      iqtree2 -s input.fasta -T $CPUS $IQTREE_EXTRA"
fi
echo " Directorio salida:      $OUTDIR/"
echo ""
echo " Reserva en Slurm:"
echo "   - Cómputo:  $CPUS CPUs (= hilos de IQ-TREE) | Memoria RAM: $MEM"
echo "   - Tiempo:   máximo $TIME | Cola: $PARTITION"
echo ""
echo " Estimación orientativa: $EST_TXT (según nº de secuencias, búsqueda, modelo y soporte)"
linea

if [ "$EST_CODE" -eq 3 ]; then
    echo ""
    echo " [AVISO] Con estas opciones el análisis probablemente no termine en el límite de 5 días."
    if ! ask_yn "¿Seguro que quieres enviarlo así?" "n"; then
        echo "Operación cancelada. Vuelve a lanzar el script y elige -fast, un modelo fijo o SH-aLRT."
        exit 0
    fi
fi

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
    --job-name="IQT_${NUM_SEQS}" \
    "$MASTER_SCRIPT" "$INPUT_FASTA" "$IQTREE_EXTRA" "$OUTDIR"

echo ""
echo " Job enviado con éxito a Slurm. ID del trabajo: $JOB_ID"
echo ""
echo " Seguimiento y gestión:"
echo "   - Estado de la cola:            squeue -u \$USER"
echo "   - Ver salida en tiempo real:     tail -f $OUTDIR/slurm_${JOB_ID}.err"
echo "   - Directorio de resultados:     $OUTDIR/"