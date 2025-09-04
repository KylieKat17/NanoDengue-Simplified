#!/bin/bash
###############################################################################
# nanoDengue (Version B – global-config edition)
# ---------------------------------------------------------------------------
# Uses one central directory for tool locations. Point the script to it with
# the env-var BIO_CONF_DIR (defaults to "$HOME/bio-configs")
#
# Expected conf file:  $BIO_CONF_DIR/nanoDengue.conf
#   Keys: MINIMAP2, SAMTOOLS, NANOPLOT, GZIP, GUNZIP
#
# All pipeline logic below this header is untouched from the OG
# two-stage workflow (NanoPlot QC → minimap2 + samtools consensus)
###############################################################################

set -euo pipefail  # -e exit on error | -u treat unset vars as error | pipefail = stricter piping

###############################################################################
# 0) LOAD PER‑MACHINE CONFIG & ALIAS BINARIES                                 #
###############################################################################

# Where to look for the global *.conf* — user can override with BIO_CONF_DIR.
CONF_ROOT="${BIO_CONF_DIR:-$HOME/bio-configs}"
CONF_FILE="$CONF_ROOT/nanoDengue.conf"

if [[ -f "$CONF_FILE" ]]; then
  echo "[INFO] Loading global tool paths from $CONF_FILE"
  # shellcheck source=/dev/null : skip static analysis (user path)
  source "$CONF_FILE"
  # Create lightweight aliases *only if* a path is set — otherwise fall back
  # to whatever is already on $PATH (should be set-up correctly in the bashrc already)
  [[ -n "${MINIMAP2:-}" ]] && alias minimap2="$MINIMAP2"
  [[ -n "${SAMTOOLS:-}"  ]] && alias samtools="$SAMTOOLS"
  [[ -n "${NANOPLOT:-}"  ]] && alias NanoPlot="$NANOPLOT"
  [[ -n "${GZIP:-}"      ]] && alias gzip="$GZIP"
  [[ -n "${GUNZIP:-}"    ]] && alias gunzip="$GUNZIP"
else
  echo "[WARN] No $CONF_FILE found. Assuming all tools are on \$PATH."
fi

###############################################################################
# 1) DEPENDENCY CHECK — VERIFY ALL REQUIRED COMMANDS EXIST                    #
# ----------------------------------------------------------------------------#
# Tests *after* aliases so custom paths count.  If something is missing,      #
# exit early with a clear message instead of failing halfway through.         #
###############################################################################
required_cmds=(minimap2 samtools NanoPlot gzip gunzip)
missing=()
for cmd in "${required_cmds[@]}"; do
  command -v "$cmd" >/dev/null 2>&1 || missing+=("$cmd")
done
if (( ${#missing[@]} )); then
  echo "[ERROR] Missing required command(s): ${missing[*]}" >&2
  echo "        Install them, add to \$PATH, or define them in $CONF_FILE" >&2
  exit 1
fi

###############################################################################
# 2) AUTO‑DETECT DEFAULT INPUT / OUTPUT PATHS                                 #
###############################################################################
# Reference FASTA: first *.fasta / *.fa in CWD, else fall back to NC_001477.  #
auto_ref=$(ls *.fasta *.fa 2>/dev/null | head -n1 || true)
DEFAULT_REF=${auto_ref:-NC_001477.fasta} 
# ^^ change if you want to, but this is the default from the OG

# Reads directory: use existing fastq_pass/ if present, else first dir containing
# any *.fastq* file.  You can *always* override with -i
if [[ -d fastq_pass ]]; then
  DEFAULT_FASTQ_DIR="fastq_pass"
else
  auto_fastq_dir=$(find . -type f -name "*.fastq*" -printf '%h\n' \
                   | sort -u | head -n1 || true)
  DEFAULT_FASTQ_DIR=${auto_fastq_dir:-fastq_pass}
fi
DEFAULT_OUT_BASE="Nano"  # where NanoPlot PNGs/HTMLs will go

# Initialise variables so flags can overwrite later
REF_GENOME="$DEFAULT_REF"
FASTQ_DIR="$DEFAULT_FASTQ_DIR"
OUT_BASE="$DEFAULT_OUT_BASE"

###############################################################################
# 3) ARG PARSING  (‑r / ‑i / ‑o)                                             #
###############################################################################
usage() {
  cat <<EOF
Usage: $0 [-r reference.fasta] [-i fastq_dir] [-o output_base] [-h]

Global config dir : $CONF_ROOT
Expected file     : nanoDengue.conf  (see header for keys)

Auto-detected defaults:
  reference.fasta : $DEFAULT_REF
  fastq_dir       : $DEFAULT_FASTQ_DIR
  output_base     : $DEFAULT_OUT_BASE

Flags override the auto-defaults (data paths), but **NOT** the global tool paths.
Tool paths still come from \$PATH or $CONF_FILE. 
(To switch tool builds, point BIO_CONF_DIR to a different directory)
EOF
}

while getopts ":r:i:o:h" opt; do
  case $opt in
    r) REF_GENOME="$OPTARG" ;;
    i) FASTQ_DIR="$OPTARG" ;;
    o) OUT_BASE="$OPTARG" ;;
    h) usage; exit 0 ;;
    \?) echo "[ERROR] Unknown option -$OPTARG" >&2; usage; exit 1 ;;
  esac
done

printf "\n  Reference genome : %s\n" "$REF_GENOME"
printf "  FASTQ directory  : %s\n" "$FASTQ_DIR"
printf "  Output base dir  : %s\n\n" "$OUT_BASE"

###############################################################################
# 4) CREATE SYMLINKS SO ORIGINAL LOOPS DON’T CARE ABOUT CUSTOM PATHS          #
# ----------------------------------------------------------------------------#
# Keeps the *original* fastq_pass/ and Nano/ inside the run folder because    #
# many downstream tools assume these names.  Symlinks are cheap & reversible. #
###############################################################################
if [[ "$FASTQ_DIR" != "fastq_pass" ]]; then
  [[ -L fastq_pass ]] && rm fastq_pass        # remove old link if exists
  ln -s "$FASTQ_DIR" fastq_pass               # create new symlink -> real dir
fi
if [[ "$OUT_BASE" != "Nano" ]]; then
  [[ -L Nano ]] && rm Nano
  ln -s "$OUT_BASE" Nano
fi
# Pass to legacy variable expected later
ref_genome="$REF_GENOME"

###############################################################################
# 5) ORIGINAL PIPELINE  (NanoPlot → minimap2 + samtools)                      #
###############################################################################

# QC step – NanoPlot
for file in fastq_pass/*/*.fastq.gz; do
  dir=$(dirname "$file")
  base=$(basename "$file" .fastq.gz)
  outdir="Nano/$dir/$base"              # QC output mirrors read subfolders
  mkdir -p "$outdir"
  NanoPlot --fastq "$file" --plots kde hex dot --outdir "$outdir"
done

# Decompress reads (makes a decompressed copy alongside original .gz)
for gz in fastq_pass/*/*.fastq.gz; do
  gunzip -c "$gz" > "${gz%.gz}"          # keep original .gz intact
done

# Alignment and consensus – minimap2 + samtools
for fq in fastq_pass/*/*.fastq; do
  sam="${fq%.fastq}_minimap2.sam"
  minimap2 -a "$ref_genome" "$fq" > "$sam"

  bam="${sam%.sam}_sorted.bam"
  samtools sort "$sam" -o "$bam"
  samtools index "$bam"

  consensus="${bam%.bam}.consensus.fasta"
  samtools consensus -f fasta "$bam" -o "$consensus"
done

###############################################################################
# End of file
###############################################################################
