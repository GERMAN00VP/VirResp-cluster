#!/bin/bash
# ==============================================================================
# FRONTAL INTERACTIVO PARA IQ-TREE 2 (Máxima Verosimilitud)
# Uso: ./lanzar_iqtree.sh alineamiento.fasta
# ==============================================================================

set -euo pipefail

# Mensaje pedagógico inicial
echo "=========================================================="
echo " LANZADOR INTERACTIVO DE IQ-TREE 2 (Clúster HPC)"
echo "=========================================================="
echo " Requisito: Debes haber alineado previamente tus secuencias (ej. con MAFFT)."
echo " Consejo: Ejecuta este script preferiblemente desde la carpeta donde está el FASTA."
echo " Los resultados se guardarán en una carpeta *_iqtree_results/"
echo "=========================================================="
echo ""

# 1. VALIDACIÓN DE ENTRADA Y RESOLUCIÓN DE RUTAS ABSOLUTAS
if [ "$#" -ne 1 ]; then
    echo " ERROR: Debes proporcionar un archivo FASTA alineado."
    echo "Uso: $0 <alineamiento.fasta>"
    exit 1
fi

if [ ! -f "$1" ]; then
    echo " ERROR: El archivo '$1' no existe."
    exit 1
fi

INPUT_FASTA="$(realpath "$1")"
WORKDIR="$(dirname "$INPUT_FASTA")"
MASTER_SCRIPT="/data/cnm/vrg/VirResp-cluster/scripts/master_iqtree.sh"

# Verificación de existencia del script ejecutor
if [ ! -x "$MASTER_SCRIPT" ]; then
    echo " ERROR CRÍTICO: No se encuentra o no es ejecutable el script maestro:"
    echo "   $MASTER_SCRIPT"
    exit 1
fi

# Sanitizar saltos de línea (\r)
sed -i 's/\r$//' "$INPUT_FASTA"

# 2. MÉTRICAS DEL ALINEAMIENTO
NUM_SEQS=$(grep -c '^>' "$INPUT_FASTA" || true)
if [ "$NUM_SEQS" -eq 0 ]; then
    echo " ERROR: El archivo no contiene cabeceras FASTA válidas ('>')."
    exit 1
fi

# Cálculo preciso de la longitud del alineamiento (incluyendo caracteres de gap/N)
ALN_LENGTH=$(awk '/^>/{next} {gsub(/[ \t\r\n]/, ""); total+=length($0)} END{if (NUM_SEQS>0) print int(total/NUM_SEQS); else print 0}' NUM_SEQS="$NUM_SEQS" "$INPUT_FASTA")

echo " ANALIZANDO ALINEAMIENTO:"
echo " -> Ruta absoluta: $INPUT_FASTA"
echo " -> Secuencias detectadas: $NUM_SEQS"
echo " -> Longitud del alineamiento: $ALN_LENGTH bp"
echo "=========================================================="

# ------------------------------------------------------------------------------
# 3. MENÚ INTERACTIVO DE MODELOS DE SUSTITUCIÓN
# ------------------------------------------------------------------------------
echo ""
echo "Selecciona el Modelo de Sustitución Evolutiva:"
echo "  1) Búsqueda automática con ModelFinder (-m MF) [Recomendado por defecto]"
echo "  2) GTR+F+I+G4   - Modelo nucleotídico completo y robusto (ADN/ARN)"
echo "  3) GTR+G        - Modelo estándar de Máxima Verosimilitud (ADN/ARN)"
echo "  4) HKY85+G      - Modelo sencillo para secuencias poco divergentes (ADN/ARN)"
echo "  5) TN93+G       - Modelo de Tamura-Nei (ADN/ARN)"
echo "  6) LG+F+G4      - Modelo estándar para Aminoácidos / Proteínas"
echo "  7) Personalizado - Escribir manualmente la sintaxis exacta de IQ-TREE"
read -p "Opción [1-7] (por defecto 1): " MODEL_OPT
MODEL_OPT=${MODEL_OPT:-1}

MODEL_SELECTED=""
case $MODEL_OPT in
    1) MODEL_SELECTED="MF" ;;
    2) MODEL_SELECTED="GTR+F+I+G4" ;;
    3) MODEL_SELECTED="GTR+G" ;;
    4) MODEL_SELECTED="HKY85+G" ;;
    5) MODEL_SELECTED="TN93+G" ;;
    6) MODEL_SELECTED="LG+F+G4" ;;
    7) 
       read -p "Introduce el nombre exacto del modelo (ej. GTR+F+R4): " CUSTOM_MODEL
       MODEL_SELECTED="$CUSTOM_MODEL"
       ;;
    *) MODEL_SELECTED="MF" ;;
esac

# ------------------------------------------------------------------------------
# 4. MENÚ INTERACTIVO DE SOPORTE DE RAMAS (BOOTSTRAP)
# ------------------------------------------------------------------------------
echo ""
echo "Método de Soporte de Ramas (Bootstrap):"
echo "  1) Ultrafast Bootstrap (UFBoot) - 1000 réplicas [Estándar y Rápido]"
echo "  2) SH-aLRT test + UFBoot        - Doble validación (1000 réplicas cada una)"
echo "  3) Sin Bootstrap                - Generar solo la topología del árbol rápido"
read -p "Opción [1-3] (por defecto 1): " BOOT_OPT
BOOT_OPT=${BOOT_OPT:-1}

BOOT_FLAGS=""
case $BOOT_OPT in
    1) BOOT_FLAGS="-bb 1000" ;;
    2) 
       BOOT_FLAGS="-alrt 1000 -bb 1000"
       if [ "$NUM_SEQS" -gt 5000 ]; then
           echo " AVISO: SH-aLRT + UFBoot con $NUM_SEQS secuencias puede demorarse varios días."
       fi
       ;;
    3) BOOT_FLAGS="" ;;
    *) BOOT_FLAGS="-bb 1000" ;;
esac

# ------------------------------------------------------------------------------
# 5. CONSTRUCCIÓN DE COMANDOS Y RECURSOS SLURM
# ------------------------------------------------------------------------------
IQTREE_EXTRA="-m $MODEL_SELECTED $BOOT_FLAGS"

PARTITION="short_idx"
TIME="12:00:00"
CPUS=16
MEM="64G"

if [ "$MODEL_SELECTED" = "MF" ] || [ "$NUM_SEQS" -gt 1000 ]; then
    PARTITION="middle_idx"
    TIME="48:00:00"
    CPUS=32
    MEM="128G"
fi

# Ajuste automático de seguridad para datasets masivos y parseo seguro de -mem
if [ "$NUM_SEQS" -gt 5000 ]; then
    echo ""
    echo " DATASET MASIVO DETECTADO ($NUM_SEQS secuencias)."
    echo "   Asignando nodo de larga duración 'long_idx' y 300 GB de RAM."
    PARTITION="long_idx"
    TIME="5-00:00:00"
    CPUS=32
    MEM="300G"
    
    # Extraer solo el número entero para pasárselo correctamente a IQ-TREE 2
    MEM_GB="${MEM%G}"
    IQTREE_EXTRA="$IQTREE_EXTRA -mem ${MEM_GB}G"
fi

echo ""
echo "=========================================================="
echo " RESUMEN DE CONFIGURACIÓN"
echo "=========================================================="
echo " -> Modelo elegido:  $MODEL_SELECTED"
echo " -> Flags IQ-TREE:   $IQTREE_EXTRA"
echo " -> Directorio:      $WORKDIR"
echo " -> Recursos Slurm:  Partición $PARTITION | $CPUS cpus | $MEM RAM | Max $TIME"
echo "=========================================================="

read -p "¿Enviar trabajo a la cola de Slurm? (s/n): " CONFIRM
if [[ "$CONFIRM" != "s" && "$CONFIRM" != "S" ]]; then
    echo "Operación cancelada por el usuario."
    exit 0
fi

# Envío a Slurm con redirección de directorio de trabajo y logs
sbatch \
  --chdir="$WORKDIR" \
  --output="$WORKDIR/iqtree_%j.out" \
  --error="$WORKDIR/iqtree_%j.err" \
  --partition="$PARTITION" \
  --time="$TIME" \
  --cpus-per-task="$CPUS" \
  --mem="$MEM" \
  --job-name="IQT_${NUM_SEQS}" \
  --export=ALL,INPUT_FASTA="$INPUT_FASTA",IQTREE_EXTRA="$IQTREE_EXTRA" \
  "$MASTER_SCRIPT"

echo ""
echo " Trabajos enviado a Slurm con éxito."
echo " Logs de seguimiento guardados en:"
echo "   $WORKDIR/iqtree_<JOB_ID>.out"
echo "   $WORKDIR/iqtree_<JOB_ID>.err"
echo " Puedes consultar el estado con: squeue -u $USER"