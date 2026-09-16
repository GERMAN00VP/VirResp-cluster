#!/bin/bash
# ==============================================================================
# FRONTAL INTERACTIVO PARA MAFFT (Alineamiento de Secuencias Masivas)
# Uso: ./lanzar_mafft.sh mi_fichero.fasta
# ==============================================================================

set -euo pipefail

# Mensaje pedagógico inicial
echo "=========================================================="
echo " LANZADOR INTERACTIVO DE MAFFT (Clúster HPC)"
echo "=========================================================="
echo " Requisito: Debes proporcionar tus secuencias en formato FASTA."
echo " Consejo: Ejecuta este script preferiblemente desde la carpeta donde está el FASTA."
echo " Los resultados se guardarán como <nombre>_aligned.fasta"
echo "=========================================================="
echo ""

# 1. VALIDACIÓN DE ENTRADA Y RESOLUCIÓN DE RUTAS ABSOLUTAS
if [ "$#" -ne 1 ]; then
    echo " ERROR: Debes proporcionar un archivo FASTA."
    echo "Uso: $0 <archivo_entrada.fasta>"
    exit 1
fi

if [ ! -f "$1" ]; then
    echo " ERROR: El archivo '$1' no existe."
    exit 1
fi

INPUT_FASTA="$(realpath "$1")"
WORKDIR="$(dirname "$INPUT_FASTA")"
MASTER_SCRIPT="/data/cnm/vrg/VirResp-cluster/scripts/master_mafft.sh"

# Verificación de existencia del script ejecutor maestro
if [ ! -x "$MASTER_SCRIPT" ]; then
    echo " ERROR CRÍTICO: No se encuentra o no es ejecutable el script maestro:"
    echo "   $MASTER_SCRIPT"
    exit 1
fi

# Sanitizar saltos de línea (\r) de Windows
sed -i 's/\r$//' "$INPUT_FASTA"

# 2. MÉTRICAS DEL ARCHIVO DE ENTRADA
NUM_SEQS=$(grep -c '^>' "$INPUT_FASTA" || true)
if [ "$NUM_SEQS" -eq 0 ]; then
    echo " ERROR: El archivo no contiene cabeceras FASTA válidas ('>')."
    exit 1
fi

# Longitud estimada de las secuencias
ALN_LENGTH=$(awk '/^>/{next} {gsub(/[ \t\r\n]/, ""); total+=length($0)} END{if (NUM_SEQS>0) print int(total/NUM_SEQS); else print 0}' NUM_SEQS="$NUM_SEQS" "$INPUT_FASTA")

echo " ANALIZANDO ARCHIVO FASTA:"
echo " -> Ruta absoluta: $INPUT_FASTA"
echo " -> Muestras/Secuencias: $NUM_SEQS"
echo " -> Longitud media estimada: $ALN_LENGTH bp"
echo "=========================================================="

# ------------------------------------------------------------------------------
# 3. CONTROLES INTERACTIVOS
# ------------------------------------------------------------------------------
echo ""
echo "Estrategias de Alineamiento (como en la Web de MAFFT):"
echo "  1) AUTO        - Recomendado. Elige el mejor algoritmo automáticamente."
echo "  2) L-INS-i     - Alta precisión (Iterativo). Ideal para < 200 seqs y regiones variables."
echo "  3) FFT-NS-2    - Rápido (Defecto estándar). Ideal para miles de secuencias."
echo "  4) PartTree    - ULTRA-MASIVO. Obligatorio si tienes > 10.000 genomas completos."
read -p "Selecciona una estrategia [1-4] (por defecto 1): " STRAT_OPT
STRAT_OPT=${STRAT_OPT:-1}

echo ""
read -p "¿Ajustar orientación de hebras en secuencias invertidas? (--adjustdirection) [S/n]: " REV_OPT
REV_OPT=${REV_OPT:-S}

echo ""
read -p "¿Reordenar secuencias alineadas según similitud? (--reorder) [s/N]: " REORD_OPT
REORD_OPT=${REORD_OPT:-N}

echo ""
read -p "¿Modificar penalización de apertura de GAPs (--op)? [Enter para valor estándar 1.53] (Aplica a ciertos algoritmos): " OP_OPT
OP_OPT=${OP_OPT:-1.53}

# ------------------------------------------------------------------------------
# 4. TRADUCCIÓN A COMANDOS Y RECURSOS SLURM
# ------------------------------------------------------------------------------
MAFFT_EXTRA=""
PARTITION="short_idx"
TIME="12:00:00"
CPUS=16
MEM="64G"

case $STRAT_OPT in
    2)
        if [ "$NUM_SEQS" -gt 500 ]; then
             echo " AVISO: L-INS-i es O(N^2) y resultará muy lento para más de 500 secuencias."
        fi
        MAFFT_EXTRA="--maxiterate 1000 --localpair"
        PARTITION="middle_idx"
        TIME="24:00:00"
        STRAT_NAME="L-INS-i (Alta Precisión)"
        ;;
    3)
        MAFFT_EXTRA="--retree 2"
        STRAT_NAME="FFT-NS-2 (Rápido)"
        ;;
    4)
        MAFFT_EXTRA="--parttree --retree 1"
        PARTITION="long_idx"
        TIME="5-00:00:00"
        CPUS=32
        MEM="300G"
        STRAT_NAME="PartTree (Ultra-masivo)"
        ;;
    *)
        MAFFT_EXTRA="--auto"
        STRAT_NAME="AUTO (Automático)"
        ;;
esac

# Modificadores adicionales
if [[ "$REV_OPT" =~ ^[Ss]$ ]]; then
    MAFFT_EXTRA="$MAFFT_EXTRA --adjustdirection"
fi

if [[ "$REORD_OPT" =~ ^[Ss]$ ]]; then
    MAFFT_EXTRA="$MAFFT_EXTRA --reorder"
fi

if [ "$OP_OPT" != "1.53" ]; then
    MAFFT_EXTRA="$MAFFT_EXTRA --op $OP_OPT"
fi

# Advertencia por combinación inusual con PartTree
if [ "$STRAT_OPT" -eq 4 ] && { [[ "$REV_OPT" =~ ^[Ss]$ ]] || [[ "$REORD_OPT" =~ ^[Ss]$ ]]; }; then
    echo " AVISO: Usar PartTree con --adjustdirection o --reorder puede provocar comportamientos anómalos."
fi

# Redirección automática para datasets masivos
if [ "$NUM_SEQS" -gt 10000 ] && [ "$STRAT_OPT" -ne 4 ]; then
    echo " Dataset masivo detectado (>10.000 seqs). Redirigiendo a nodo grande en 'long_idx'."
    PARTITION="long_idx"
    TIME="5-00:00:00"
    CPUS=32
    MEM="300G"
fi

echo ""
echo "=========================================================="
echo " RESUMEN DEL TRABAJO A ENVIAR"
echo "=========================================================="
echo " -> Estrategia: $STRAT_NAME"
echo " -> Flags de MAFFT: $MAFFT_EXTRA"
echo " -> Directorio: $WORKDIR"
echo " -> Recursos Slurm: Partición $PARTITION | $CPUS cpus | $MEM RAM | Max $TIME"
echo "=========================================================="

read -p "¿Lanzar al clúster? (s/n): " CONFIRM
if [[ "$CONFIRM" != "s" && "$CONFIRM" != "S" ]]; then
    echo "Operación cancelada por el usuario."
    exit 0
fi

# Envío a Slurm fijando la ruta del usuario para la salida de logs
sbatch \
  --chdir="$WORKDIR" \
  --output="$WORKDIR/mafft_%j.out" \
  --error="$WORKDIR/mafft_%j.err" \
  --partition="$PARTITION" \
  --time="$TIME" \
  --cpus-per-task="$CPUS" \
  --mem="$MEM" \
  --job-name="MAFFT_${NUM_SEQS}" \
  --export=ALL,INPUT_FASTA="$INPUT_FASTA",MAFFT_EXTRA="$MAFFT_EXTRA",NUM_SEQS="$NUM_SEQS" \
  "$MASTER_SCRIPT"

echo ""
echo " Trabajo enviado con éxito."
echo " Logs de seguimiento guardados en:"
echo "   $WORKDIR/mafft_<JOB_ID>.out"
echo "   $WORKDIR/mafft_<JOB_ID>.err"
echo " Consulta el estado con: squeue -u $USER"