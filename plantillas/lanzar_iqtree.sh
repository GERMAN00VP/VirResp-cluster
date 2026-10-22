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

if [ "$NUM_SEQS" -gt 200 ]; then DEF_THR=2; else DEF_THR=1; fi
echo ""
echo "3. HILOS DE CÓMPUTO"
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
echo "4. FIJACIÓN DE SEMILLA ALEATORIA PARA REPRODUCIBILIDAD (--seed)"
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

# Con AUTO el lanzador pasa -T y el maestro no lo añade
if [ -n "$THREAD_FLAGS" ]; then
    FLAGS+=($THREAD_FLAGS)
fi

# Factor de esfuerzo para Slurm
WORK=$(( NUM_SEQS * ALN_LEN ))
IS_MF=0
if [[ "$MODEL_FLAG" == *"-m MF"* ]]; then IS_MF=1; fi

if [ "$NUM_SEQS" -gt 5000 ] || [ "$WORK" -gt 2000000000 ]; then
    PARTITION="long_idx";   TIME="5-00:00:00"; MEM_GB=300
elif [ "$NUM_SEQS" -gt 1000 ] || [ "$IS_MF" -eq 1 ] || [ "$WORK" -gt 100000000 ]; then
    PARTITION="middle_idx"; TIME="48:00:00";   MEM_GB=128
else
    PARTITION="short_idx";  TIME="12:00:00";   MEM_GB=64
fi

# Añadir límite de memoria explícito a IQ-TREE 2
FLAGS+=("-mem" "${MEM_GB}G")
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
    --job-name="IQT_${NUM_SEQS}" \
    "$MASTER_SCRIPT" "$INPUT_FASTA" "$IQTREE_EXTRA" "$OUTDIR"

echo ""
echo " Job enviado con éxito a Slurm. ID del trabajo: $JOB_ID"
echo ""
echo " Seguimiento y gestión:"
echo "   - Estado de la cola:            squeue -u \$USER"
echo "   - Ver salida en tiempo real:     tail -f $OUTDIR/slurm_${JOB_ID}.err"
echo "   - Directorio de resultados:     $OUTDIR/"