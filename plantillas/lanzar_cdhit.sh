#!/bin/bash
# ==============================================================================
# FRONTAL INTERACTIVO PARA CD-HIT (Desduplicación y Haplotipado de Virus)
# Uso: ./lanzar_cdhit.sh mi_fichero.fasta
# ==============================================================================

set -euo pipefail

# Mensaje pedagógico inicial
echo "=========================================================="
echo " LANZADOR INTERACTIVO DE CD-HIT (Clúster HPC)"
echo "=========================================================="
echo " Requisito: Debes proporcionar un archivo FASTA."
echo " Uso principal: Colapsar secuencias 100% idénticas (haplotipos)."
echo " Nota: Funciona de forma óptima con secuencias sin alinear."
echo " Resultados: FASTA desduplicado y reporte .haplotypes.txt"
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
MASTER_SCRIPT="/data/cnm/vrg/VirResp-cluster/scripts/master_cdhit.sh"

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

ALN_LENGTH=$(awk '/^>/{next} {gsub(/[ \t\r\n]/, ""); total+=length($0)} END{if (NUM_SEQS>0) print int(total/NUM_SEQS); else print 0}' NUM_SEQS="$NUM_SEQS" "$INPUT_FASTA")

echo " ANALIZANDO ARCHIVO FASTA:"
echo " -> Ruta absoluta: $INPUT_FASTA"
echo " -> Muestras/Secuencias: $NUM_SEQS"
echo " -> Longitud media estimada: $ALN_LENGTH bp"
echo "=========================================================="

# ------------------------------------------------------------------------------
# 3. CONTROLES INTERACTIVOS (SELECCIÓN DE SIMILITUD Y PARÁMETROS)
# ------------------------------------------------------------------------------
echo ""
echo "Selecciona el Umbral de Similitud:"
echo "  1) 100% Identidad (-c 1.0)  - Desduplicación pura (Haplotipado estricto)"
echo "  2)  99% Identidad (-c 0.99) - Colapsar secuencias con errores de lectura o microvariación"
echo "  3)  95% Identidad (-c 0.95) - Agrupar por linajes/subvariantes genómicas"
echo "  4) Personalizado            - Umbral específico (ej. 0.90)"
read -p "Opción [1-4] (por defecto 1): " SIM_OPT
SIM_OPT=${SIM_OPT:-1}

ID_THRESHOLD="1.0"
case $SIM_OPT in
    1) ID_THRESHOLD="1.0" ;;
    2) ID_THRESHOLD="0.99" ;;
    3) ID_THRESHOLD="0.95" ;;
    4) 
       read -p "Introduce el umbral de identidad (0.80 - 1.0): " CUSTOM_ID
       # Validación matemática del rango con awk
       if ! awk "BEGIN{exit !($CUSTOM_ID >= 0.80 && $CUSTOM_ID <= 1.0)}"; then
           echo " ERROR: El umbral personalizado debe estar entre 0.80 y 1.0"
           exit 1
       fi
       ID_THRESHOLD="$CUSTOM_ID"
       ;;
    *) ID_THRESHOLD="1.0" ;;
esac

# Asignación del tamaño del word (-n) requerida técnicamente por CD-HIT
WORD_SIZE=10
if awk "BEGIN{exit !($ID_THRESHOLD < 0.88)}"; then
    WORD_SIZE=5
elif awk "BEGIN{exit !($ID_THRESHOLD < 0.90)}"; then
    WORD_SIZE=6
elif awk "BEGIN{exit !($ID_THRESHOLD < 0.92)}"; then
    WORD_SIZE=8
else
    WORD_SIZE=10
fi

echo ""
read -p "Cobertura mínima de la secuencia más corta (-aS) [Enter para 1.0 (100%)]: " AS_OPT
AS_OPT=${AS_OPT:-1.0}

# -d 0 previene que CD-HIT corte los IDs largos de las secuencias a 20 caracteres
CDHIT_EXTRA="-c $ID_THRESHOLD -n $WORD_SIZE -aS $AS_OPT -g 1 -d 0"

# ------------------------------------------------------------------------------
# 4. CONFIGURACIÓN DE RECURSOS SLURM
# ------------------------------------------------------------------------------
PARTITION="short_idx"
TIME="04:00:00"
CPUS=8
MEM="32G"

if [ "$NUM_SEQS" -gt 50000 ]; then
    echo " Dataset masivo (>50.000 seqs). Escalando asignación en cluster."
    PARTITION="middle_idx"
    TIME="24:00:00"
    CPUS=16
    MEM="64G"
fi

echo ""
echo "=========================================================="
echo " RESUMEN DEL TRABAJO A ENVIAR"
echo "=========================================================="
echo " -> Objetivo:          $( [ "$ID_THRESHOLD" = "1.0" ] && echo "Haplotipado / Desduplicación estricta" || echo "Agrupamiento al $ID_THRESHOLD" )"
echo " -> Flags CD-HIT:      $CDHIT_EXTRA"
echo " -> Directorio:        $WORKDIR"
echo " -> Recursos Slurm:    Partición $PARTITION | $CPUS cpus | $MEM RAM | Max $TIME"
echo "=========================================================="

read -p "¿Lanzar al clúster? (s/n): " CONFIRM
if [[ "$CONFIRM" != "s" && "$CONFIRM" != "S" ]]; then
    echo "Operación cancelada por el usuario."
    exit 0
fi

sbatch \
  --chdir="$WORKDIR" \
  --output="$WORKDIR/cdhit_%j.out" \
  --error="$WORKDIR/cdhit_%j.err" \
  --partition="$PARTITION" \
  --time="$TIME" \
  --cpus-per-task="$CPUS" \
  --mem="$MEM" \
  --job-name="CDHIT_${NUM_SEQS}" \
  --export=ALL,INPUT_FASTA="$INPUT_FASTA",CDHIT_EXTRA="$CDHIT_EXTRA" \
  "$MASTER_SCRIPT"

echo ""
echo " Trabajo enviado con éxito."
echo " Logs de ejecución guardados en:"
echo "   $WORKDIR/cdhit_<JOB_ID>.out"
echo "   $WORKDIR/cdhit_<JOB_ID>.err"
echo " Consulta el estado con: squeue -u $USER"