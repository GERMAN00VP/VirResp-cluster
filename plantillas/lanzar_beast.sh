#!/bin/bash
# ==============================================================================
# FRONTAL INTERACTIVO PARA BEAST + BEAGLE (Filodinámica y Relojes Moleculares)
# Uso: ./lanzar_beast.sh mi_analisis.xml
# ==============================================================================

set -euo pipefail

# Mensaje pedagógico inicial para orientar al virólogo
echo "=========================================================="
echo "LANZADOR INTERACTIVO DE BEAST + BEAGLE (Clúster HPC)"
echo "=========================================================="
echo " Requisito: Debes haber generado tu archivo .xml previamente en BEAUti."
echo " Consejo: Ejecuta este script desde la carpeta donde tienes tu XML."
echo " Los resultados se guardarán en una carpeta *_beast_results/"
echo "=========================================================="
echo ""

# Validar argumento de entrada
if [ "$#" -ne 1 ]; then
    echo "ERROR: Debes proporcionar un archivo XML."
    echo "Uso: $0 <analisis.xml>"
    exit 1
fi

# 1. RESOLUCIÓN DE RUTAS ABSOLUTAS (Fix 1 y Fix 3)
if [ ! -f "$1" ]; then
    echo "ERROR: El archivo '$1' no existe."
    exit 1
fi

INPUT_XML="$(realpath "$1")"
WORKDIR="$(dirname "$INPUT_XML")"
MASTER_SCRIPT="/data/cnm/vrg/scripts/master_beast.sh"

# Verificar existencia y permisos del script maestro ejecutor (Fix 8)
if [ ! -x "$MASTER_SCRIPT" ]; then
    echo "ERROR CRÍTICO: No se encuentra o no es ejecutable el script maestro:"
    echo "   $MASTER_SCRIPT"
    echo "   Asegúrate de que la ruta sea correcta y tenga permisos (+x)."
    exit 1
fi

# Sanitizar saltos de línea (\r) de Windows
sed -i 's/\r$//' "$INPUT_XML"

# 2. PARSEO ROBUSTO DEL CHAINLENGTH (Fix 2)
CHAIN_LENGTH=$(grep -m1 -oP 'mcmc[^>]*chainLength="\K[0-9]+' "$INPUT_XML" || echo "Desconocida")

echo "ARCHIVO XML IDENTIFICADO:"
echo " -> Ruta absoluta: $INPUT_XML"
echo " -> Longitud de cadena MCMC (chainLength): $CHAIN_LENGTH"
echo "=========================================================="

# ------------------------------------------------------------------------------
# MENÚ DE SELECCIÓN DE MOTOR COMPUTACIONAL
# ------------------------------------------------------------------------------
echo ""
echo "Selecciona el Motor de Aceleración Computacional:"
echo "  1) CPU Multinúcleo + BEAGLE (Estándar y compatible) [Recomendado]"
echo "  2) GPU + BEAGLE (Alta velocidad para matrices masivas o genomas completos)"
echo "  3) BEAST Nativo (Sin aceleración BEAGLE)"
read -p "Opción [1-3] (por defecto 1): " ENGINE_OPT
ENGINE_OPT=${ENGINE_OPT:-1}

BEAGLE_FLAGS=""
USE_GPU="N"

case $ENGINE_OPT in
    1)
        BEAGLE_FLAGS="-beagle -beagle_CPU -beagle_SSE"
        ENGINE_NAME="BEAGLE (CPU/SSE)"
        ;;
    2)
        BEAGLE_FLAGS="-beagle -beagle_GPU"
        ENGINE_NAME="BEAGLE (GPU)"
        USE_GPU="Y"
        ;;
    3)
        BEAGLE_FLAGS=""
        ENGINE_NAME="BEAST Nativo"
        ;;
    *)
        BEAGLE_FLAGS="-beagle -beagle_CPU -beagle_SSE"
        ENGINE_NAME="BEAGLE (CPU/SSE)"
        ;;
esac

# ------------------------------------------------------------------------------
# PRECISIÓN NUMÉRICA CON PROTECCIÓN PARA GPU (Fix 9)
# ------------------------------------------------------------------------------
PRECISION_FLAG=""

if [ "$USE_GPU" = "Y" ]; then
    # En la mayoría de GPUs comerciales/HPC la doble precisión colapsa el rendimiento o falla
    PRECISION_FLAG="-beagle_single"
    echo ""
    echo "AVISO: Aceleración GPU detectada. Se fuerza precisión simple (-beagle_single)"
    echo "   para evitar problemas de compatibilidad y pérdida drástica de velocidad."
elif [ "$ENGINE_OPT" -ne 3 ]; then
    echo ""
    echo "Selecciona la Precisión de Cálculo de BEAGLE:"
    echo "  1) Doble Precisión (-beagle_double) [Estabilidad en filodinámica/relojes complejos]"
    echo "  2) Precisión Simple (-beagle_single) [Más rápido, riesgo numérico en cadenas muy largas]"
    read -p "Opción [1-2] (por defecto 1): " PREC_OPT
    PREC_OPT=${PREC_OPT:-1}

    if [ "$PREC_OPT" -eq 2 ]; then
        PRECISION_FLAG="-beagle_single"
    else
        PRECISION_FLAG="-beagle_double"
    fi
fi

# ------------------------------------------------------------------------------
# ASIGNACIÓN DE RECURSOS SLURM
# ------------------------------------------------------------------------------
PARTITION="middle_idx"
TIME="48:00:00"
CPUS=16
MEM="64G"
GPU_SLURM_DIRECTIVE=""

if [ "$USE_GPU" = "Y" ]; then
    PARTITION="gpu"
    TIME="5-00:00:00"
    CPUS=8
    MEM="64G"
    GPU_SLURM_DIRECTIVE="--gres=gpu:1"
else
    echo ""
    echo "Estimación de tiempo de ejecución para CPU:"
    echo "  1) MCMC Corto/Medio (Hasta 48 horas)"
    echo "  2) MCMC Largo/Muy pesado (Hasta 5 días - Nodos de larga duración)"
    read -p "Opción [1-2] (por defecto 1): " TIME_OPT
    TIME_OPT=${TIME_OPT:-1}

    if [ "$TIME_OPT" -eq 2 ]; then
        PARTITION="long_idx"
        TIME="5-00:00:00"
        CPUS=32
        MEM="128G"
    fi
fi

BEAST_EXTRA="$BEAGLE_FLAGS $PRECISION_FLAG"

echo ""
echo "=========================================================="
echo "RESUMEN DEL TRABAJO A ENVIAR A SLURM"
echo "=========================================================="
echo " -> Motor/Aceleración: $ENGINE_NAME"
echo " -> Flags BEAST/BEAGLE: $BEAST_EXTRA"
echo " -> Directorio trabajo: $WORKDIR"
echo " -> Recursos Slurm:    Partición $PARTITION | $CPUS CPUs | $MEM RAM | Tiempo máx $TIME"
echo "=========================================================="

read -p "¿Lanzar análisis MCMC al clúster? (s/n): " CONFIRM
if [[ "$CONFIRM" != "s" && "$CONFIRM" != "S" ]]; then
    echo "Operación cancelada por el usuario."
    exit 0
fi

# ------------------------------------------------------------------------------
# ENVÍO A SLURM CON FIXES 3 Y 4 (--chdir y redirection de logs por trabajo)
# ------------------------------------------------------------------------------
sbatch \
  --chdir="$WORKDIR" \
  --output="$WORKDIR/beast_%j.out" \
  --error="$WORKDIR/beast_%j.err" \
  --partition="$PARTITION" \
  --time="$TIME" \
  --cpus-per-task="$CPUS" \
  --mem="$MEM" \
  $GPU_SLURM_DIRECTIVE \
  --job-name="BEAST_MCMC" \
  --export=ALL,INPUT_XML="$INPUT_XML",BEAST_EXTRA="$BEAST_EXTRA" \
  "$MASTER_SCRIPT"

echo ""
echo "✔ Trabajos enviado a Slurm con éxito."
echo " Logs de seguimiento en tiempo real guardados en:"
echo "   $WORKDIR/beast_<JOB_ID>.out"
echo "   $WORKDIR/beast_<JOB_ID>.err"
echo " Puedes monitorear la cola con: squeue -u $USER"