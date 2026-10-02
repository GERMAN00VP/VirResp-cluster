#!/bin/bash
# ==============================================================================
# FRONTAL AUTOMÁTICO PARA BEAST + BEAGLE (Filodinámica y Relojes Moleculares)
# Uso: ./lanzar_beast.sh mi_analisis.xml
# ==============================================================================

set -euo pipefail

echo "=========================================================="
echo " LANZADOR AUTOMÁTICO DE BEAST + BEAGLE (Clúster HPC)"
echo "=========================================================="
echo " Requisito: Debes haber generado tu archivo .xml en BEAUti."
echo " Los resultados se guardarán en *_beast_results/"
echo "=========================================================="
echo ""

# 1. VALIDACIÓN DE ENTRADA Y RUTAS ABSOLUTAS
if [ "$#" -ne 1 ]; then
    echo " ERROR: Debes proporcionar un archivo XML."
    echo "Uso: $0 <analisis.xml>"
    exit 1
fi

if [ ! -f "$1" ]; then
    echo " ERROR: El archivo '$1' no existe."
    exit 1
fi

INPUT_XML="$(realpath "$1")"
WORKDIR="$(dirname "$INPUT_XML")"
MASTER_SCRIPT="/data/cnm/vrg/VirResp-cluster/scripts/master_beast.sh"

if [ ! -x "$MASTER_SCRIPT" ]; then
    echo " ERROR CRÍTICO: No se encuentra o no es ejecutable el script maestro:"
    echo "   $MASTER_SCRIPT"
    exit 1
fi

# Sanitizar saltos de línea (\r) de Windows
sed -i 's/\r$//' "$INPUT_XML"

# 2. INSPECCIÓN AUTOMÁTICA DEL XML
CHAIN_LENGTH=$(grep -m1 -oP 'mcmc[^>]*chainLength="\K[0-9]+' "$INPUT_XML" || echo "0")
NUM_SEQS=$(grep -c '<sequence>' "$INPUT_XML" || grep -c '<taxon id=' "$INPUT_XML" || echo "0")

echo " ANALIZANDO ARCHIVO XML:"
echo " -> Ruta absoluta: $INPUT_XML"
echo " -> MCMC chainLength: ${CHAIN_LENGTH:-Desconocida}"
echo " -> Taxones/Secuencias estimadas: $NUM_SEQS"
echo "=========================================================="

# 3. ASIGNACIÓN INTELIGENTE DE RECURSOS SLURM
# Configuración predeterminada estable (CPU + BEAGLE SSE + Doble Precisión)
BEAST_EXTRA="-beagle -beagle_CPU -beagle_SSE -beagle_double"

PARTITION="middle_idx"
TIME="48:00:00"
CPUS=16
MEM="64G"
REASON="Análisis estándar (Cadena <= 50M)"

# Escalar a la partición de larga duración si la cadena MCMC es muy grande o hay muchas secuencias
if [ "$CHAIN_LENGTH" -gt 50000000 ] 2>/dev/null || [ "$NUM_SEQS" -gt 300 ] 2>/dev/null; then
    PARTITION="long_idx"
    TIME="5-00:00:00"
    CPUS=32
    MEM="128G"
    REASON="Dataset masivo o MCMC larga (>50M)"
fi

echo ""
echo "=========================================================="
echo " RESUMEN DE LA CONFIGURACIÓN CALCULADA"
echo "=========================================================="
echo " -> Motor de Cálculo: BEAGLE CPU (Doble Precisión)"
echo " -> Causa de asignación: $REASON"
echo " -> Directorio trabajo: $WORKDIR"
echo " -> Recursos Slurm:    Partición $PARTITION | $CPUS CPUs | $MEM RAM | Max $TIME"
echo "=========================================================="

read -p "¿Lanzar análisis MCMC al clúster? (s/n): " CONFIRM
if [[ "$CONFIRM" != "s" && "$CONFIRM" != "S" ]]; then
    echo "Operación cancelada por el usuario."
    exit 0
fi

# 4. ENVÍO A SLURM
sbatch \
  --chdir="$WORKDIR" \
  --output="$WORKDIR/beast_%j.out" \
  --error="$WORKDIR/beast_%j.err" \
  --partition="$PARTITION" \
  --time="$TIME" \
  --cpus-per-task="$CPUS" \
  --mem="$MEM" \
  --job-name="BEAST_MCMC" \
  --export=ALL,INPUT_XML="$INPUT_XML",BEAST_EXTRA="$BEAST_EXTRA" \
  "$MASTER_SCRIPT"

echo ""
echo " Trabajos enviado a Slurm con éxito."
echo " Logs de seguimiento en tiempo real guardados en:"
echo "   $WORKDIR/beast_<JOB_ID>.out"
echo "   $WORKDIR/beast_<JOB_ID>.err"
echo " Puedes monitorear la cola con: squeue -u $USER"