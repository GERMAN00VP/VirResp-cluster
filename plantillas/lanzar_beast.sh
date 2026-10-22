#!/bin/bash
# ==============================================================================
# LANZADOR AUTOMÁTICO DE BEAST + BEAGLE (VirResp-cluster)
# Uso: ./lanzar_beast.sh mi_analisis.xml
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
echo " FILODINÁMICA Y RELOJES MOLECULARES CON BEAST + BEAGLE"
linea
echo " Fichero de entrada: Configuración MCMC en formato XML (*.xml)"
echo " Nota: Toda la parametrización evolutiva debe haberse configurado en BEAUti."
linea

# ------------------------------------------------------------------------------
# 1. COMPROBACIONES PREVIAS
# ------------------------------------------------------------------------------
if [ "$#" -ne 1 ]; then
    error "Indica el archivo XML de configuración de BEAST." "Uso: $0 <analisis.xml>"
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

INPUT_XML="$(realpath "$1")"
WORKDIR="$(dirname "$INPUT_XML")"

if [ ! -w "$WORKDIR" ]; then
    error "Sin permiso de escritura en la carpeta:" "$WORKDIR"
fi

MASTER_SCRIPT="$(realpath -m "$SCRIPT_DIR/../scripts/master_beast.sh")"

if [ ! -x "$MASTER_SCRIPT" ]; then
    error "No existe o no es ejecutable el script maestro:" "$MASTER_SCRIPT"
fi

# Sanitizar saltos de línea de Windows (\r)
sed -i 's/\r$//' "$INPUT_XML"

# ------------------------------------------------------------------------------
# 2. INSPECCIÓN AUTOMÁTICA DEL ARCHIVO XML
# ------------------------------------------------------------------------------
if ! grep -q -E '<beast|<mcmc' "$INPUT_XML"; then
    error "El archivo '$INPUT_XML' no parece un archivo de configuración válido de BEAST."
fi

CHAIN_LENGTH=$(grep -m1 -oP '(?:chainLength|chain_length)="\K[0-9]+' "$INPUT_XML" || echo "0")
NUM_SEQS=$(grep -c -E '<sequence|<taxon id=' "$INPUT_XML" || echo "0")

echo " INSPECCIÓN DE CONFIGURACIÓN XML"
echo " -> Archivo:            $INPUT_XML"
if [ "$CHAIN_LENGTH" -gt 0 ]; then
    echo " -> Longitud cadena MCMC: $CHAIN_LENGTH iteraciones"
else
    echo " -> Longitud cadena MCMC: No detectada automáticamente"
fi
echo " -> Taxones/Secuencias:  $NUM_SEQS detectadas"
linea

# ------------------------------------------------------------------------------
# 3. MOTOR DE CÁLCULO Y ESTIMACIÓN DE RECURSOS SLURM
# ------------------------------------------------------------------------------
# Configuración predeterminada optimizada con BEAGLE (CPU + SSE + Doble Precisión)
BEAST_EXTRA="-beagle -beagle_CPU -beagle_SSE -beagle_double"

if [ "$CHAIN_LENGTH" -gt 50000000 ] || [ "$NUM_SEQS" -gt 300 ]; then
    PARTITION="long_idx"
    TIME="5-00:00:00"
    DEF_HILOS=32
    MEM="128G"
    REASON="Dataset masivo o MCMC muy larga (> 50M iteraciones)"
elif [ "$CHAIN_LENGTH" -le 10000000 ] && [ "$NUM_SEQS" -le 100 ] && [ "$CHAIN_LENGTH" -gt 0 ]; then
    PARTITION="short_idx"
    TIME="12:00:00"
    DEF_HILOS=16
    MEM="64G"
    REASON="Dataset pequeño / MCMC corta (<= 10M iteraciones)"
else
    PARTITION="middle_idx"
    TIME="48:00:00"
    DEF_HILOS=16
    MEM="64G"
    REASON="Dataset estándar (MCMC <= 50M iteraciones)"
fi

pedir_hilos "HILOS DE CÓMPUTO (-threads)" "$DEF_HILOS" \
    "BEAGLE reparte la verosimilitud entre hilos; con pocas particiones más hilos no acelera." \
    "Se reservarán tantas CPUs como hilos."
CPUS="$HILOS"

BASENAME="$(basename "$INPUT_XML" .xml)"
pedir_outdir "$WORKDIR/${BASENAME}_beast_results" "$WORKDIR"

# ------------------------------------------------------------------------------
# 4. RESUMEN Y CONFIRMACIÓN
# ------------------------------------------------------------------------------
echo ""
linea
echo " RESUMEN DE LA CONFIGURACIÓN DE BEAST"
linea
echo " Archivo de entrada:     $INPUT_XML"
echo " Motor de aceleración:   BEAGLE CPU (Aceleración SSE, Doble Precisión)"
echo " Criterio de recursos:   $REASON"
echo " Comand line flags:      beast -threads $CPUS $BEAST_EXTRA input.xml"
echo " Directorio salida:      $OUTDIR/"
echo ""
echo " Reserva en Slurm:"
echo "   - Cómputo:  $CPUS CPUs | Memoria RAM: $MEM"
echo "   - Tiempo:   máximo $TIME | Cola: $PARTITION"
linea

if ! ask_yn "¿Confirmar y enviar análisis MCMC a Slurm?" "s"; then
    echo "Operación cancelada por el usuario."
    exit 0
fi

# ------------------------------------------------------------------------------
# 5. ENVÍO DEL JOB A SLURM
# ------------------------------------------------------------------------------
enviar_job \
    --chdir="$WORKDIR" \
    --output="$OUTDIR/slurm_%j.out" \
    --error="$OUTDIR/slurm_%j.err" \
    --partition="$PARTITION" \
    --time="$TIME" \
    --cpus-per-task="$CPUS" \
    --mem="$MEM" \
    --job-name="BEAST_MCMC" \
    "$MASTER_SCRIPT" "$INPUT_XML" "$BEAST_EXTRA" "$OUTDIR"

echo ""
echo " Job enviado con éxito a Slurm. ID del trabajo: $JOB_ID"
echo ""
echo " Seguimiento y gestión:"
echo "   - Estado de la cola:            squeue -u \$USER"
echo "   - Ver salida en tiempo real:     tail -f $OUTDIR/slurm_${JOB_ID}.out"
echo "   - Directorio de resultados:     $OUTDIR/"