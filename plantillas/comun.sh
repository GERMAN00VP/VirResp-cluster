#!/bin/bash
# ==============================================================================
# FUNCIONES COMUNES DE LOS LANZADORES (VirResp-cluster)
# Uso: source "$SCRIPT_DIR/comun.sh"   (no se ejecuta directamente)
# ==============================================================================

# ------------------------------------------------------------------------------
# Funciones auxiliares de interfaz
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
# Comprobaciones rápidas de FASTA (por muestreo; no recorren el fichero entero)
# ------------------------------------------------------------------------------

# Número de secuencias, primer carácter y saltos de línea Windows.
# Define: NUM_SEQS, HAS_CR (1/0)
fasta_basico() {
    local f="$1"
    if [ "$(head -c 1 "$f")" != ">" ]; then
        error "El archivo '$f' no empieza por '>'." "No parece un FASTA válido."
    fi
    NUM_SEQS=$(grep -c '^>' "$f" || true)
    HAS_CR=0
    if head -c 100000 "$f" | grep -q $'\r'; then
        HAS_CR=1
        echo " [AVISO] El archivo tiene saltos de línea de Windows (\\r). Se tolerarán en el análisis."
    fi
}

# Longitud mínima y máxima de las primeras N secuencias y longitud de la última.
# Define: MIN_LEN, MAX_LEN, N_MUESTRA, LAST_LEN, AVG_LEN (media de la muestra)
fasta_longitudes_muestra() {
    local f="$1" n="${2:-50}" tail_bytes
    # awk lee el fichero directamente y sale al llegar a la secuencia N+1 (sin SIGPIPE)
    read -r MIN_LEN MAX_LEN N_MUESTRA AVG_LEN < <(awk -v N="$n" '
        /^>/ { if (n) L[n] = len; if (n == N) { done = 1; exit } n++; len = 0; next }
        { gsub(/[ \t\r]/, ""); len += length($0) }
        END {
            if (!done && n) L[n] = len
            min = max = L[1]; tot = 0
            for (i = 1; i <= n; i++) { if (L[i] < min) min = L[i]; if (L[i] > max) max = L[i]; tot += L[i] }
            print min + 0, max + 0, n + 0, (n ? int(tot / n) : 0)
        }' "$f")

    # Última secuencia: tail -c hace seek, no recorre todo el fichero
    tail_bytes=$(( MAX_LEN * 3 + 1000000 ))
    LAST_LEN=$(tail -c "$tail_bytes" "$f" | awk '
        /^>/ { len = 0; next }
        { gsub(/[ \t\r]/, ""); len += length($0) }
        END { print len + 0 }')
}

# Comprueba que la muestra de longitudes es homogénea (alineamiento).
# Define: ALN_LEN
fasta_exigir_alineado() {
    if [ "$MIN_LEN" -ne "$MAX_LEN" ] || [ "$LAST_LEN" -ne "$MAX_LEN" ]; then
        error "Longitudes desiguales en la muestra ($MIN_LEN-$MAX_LEN; última: $LAST_LEN)." \
              "El archivo no parece alineado. Ejecuta primero 'lanzar_mafft.sh'."
    fi
    ALN_LEN="$MAX_LEN"
}

# Tipo de secuencia a partir de una muestra de ~2 MB del inicio.
# Define: NT_PCT (porcentaje de ACGTUN sobre residuos no hueco)
fasta_tipo_muestra() {
    local f="$1" sample total nuc
    # head primero: evita SIGPIPE con pipefail
    sample=$(head -c 2000000 "$f" | grep -v '^>' | tr -d '\r\n\t .?-' | tr 'a-z' 'A-Z' || true)
    total=${#sample}
    [ "$total" -gt 0 ] || error "La muestra del archivo solo contiene huecos o está vacía."
    nuc=$(printf '%s' "$sample" | tr -cd 'ACGTUN' | wc -c)
    NT_PCT=$(( 100 * nuc / total ))
}

# Primeros identificadores duplicados (vacío si no hay)
fasta_ids_duplicados() {
    grep '^>' "$1" | sed 's/^>//; s/[[:space:]].*$//' | sort | uniq -d | head -n 5 || true
}

# Número de identificadores con caracteres problemáticos en Newick
fasta_ids_newick() {
    grep -c -E '^>[^[:space:]]*[][(),:;]' "$1" || true
}

# ------------------------------------------------------------------------------
# Directorio de salida elegible y versionado (_v2, _v3...)
# ------------------------------------------------------------------------------

# Uso: OUTDIR=$(resolver_outdir "<entrada_usuario>" "<valor_por_defecto>" "<carpeta_del_input>")
#   Intro          -> valor por defecto
#   /ruta/absoluta -> esa ruta
#   ~/algo         -> $HOME/algo
#   otra cosa      -> relativa a la carpeta del archivo de entrada
# Los avisos van a stderr para no contaminar la salida capturada.
resolver_outdir() {
    local entrada="$1" defecto="$2" base_input="$3" base candidato n=2
    entrada="$(echo "$entrada" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')"
    if [ -z "$entrada" ]; then
        base="$defecto"
    elif [[ "$entrada" == "~"* ]]; then
        base="$(realpath -m "${HOME}${entrada#\~}")"
    elif [[ "$entrada" == /* ]]; then
        base="$(realpath -m "$entrada")"
    else
        base="$(realpath -m "$base_input/$entrada")"
    fi
    base="${base%/}"
    candidato="$base"
    while [ -e "$candidato" ]; do
        candidato="${base}_v${n}"
        n=$(( n + 1 ))
    done
    if [ "$candidato" != "$base" ]; then
        echo " [AVISO] '$base' ya existe; se usará '$candidato'." >&2
    fi
    echo "$candidato"
}

# Pregunta el directorio de salida hasta obtener uno válido.
# Uso: pedir_outdir "<valor_por_defecto>" "<carpeta_del_input>"   Define: OUTDIR
pedir_outdir() {
    local defecto="$1" base_input="$2" entrada limpio padre
    echo ""
    echo " DIRECTORIO DE SALIDA"
    echo "   Ruta absoluta (/...) o ~/...: se usa tal cual."
    echo "   Nombre o ruta relativa: se crea junto al archivo de entrada ($base_input)."
    echo "   Si ya existe, se añade _v2, _v3... (nunca se sobrescribe)."
    while true; do
        read -r -p "Directorio de salida (Intro = $defecto): " entrada || { echo ""; echo "Entrada cancelada."; exit 1; }
        limpio="$(echo "$entrada" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')"
        if [ -n "$entrada" ] && [ -z "$limpio" ]; then
            echo "   Nombre vacío: escribe un nombre o pulsa Intro."
            continue
        fi
        if [[ "$limpio" == *[[:space:]]* ]]; then
            echo "   El nombre no puede contener espacios."
            continue
        fi
        OUTDIR="$(resolver_outdir "$limpio" "$defecto" "$base_input")"
        # Carpeta existente más cercana: debe ser un directorio con permiso de escritura
        padre="$(dirname "$OUTDIR")"
        while [ ! -e "$padre" ]; do padre="$(dirname "$padre")"; done
        if [ ! -d "$padre" ] || [ ! -w "$padre" ]; then
            echo "   Sin permiso de escritura para crear '$OUTDIR' (en '$padre')."
            continue
        fi
        return 0
    done
}

# Crea OUTDIR (reserva el nombre), envía el job y lo borra si sbatch falla.
# Uso: enviar_job <opciones de sbatch y script con argumentos>   Define: JOB_ID
enviar_job() {
    local salida
    mkdir -p "$OUTDIR" || error "No se pudo crear el directorio de salida:" "$OUTDIR"
    if ! salida=$(sbatch --parsable "$@"); then
        rmdir "$OUTDIR" 2>/dev/null || true
        error "sbatch no pudo enviar el trabajo."
    fi
    JOB_ID="${salida%%;*}"
}
