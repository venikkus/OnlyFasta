#!/usr/bin/env bash
set -euo pipefail

# Resolve the repository from this wrapper location; dataset and output paths come from config/config.yaml.
PROJECT_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
SNAKEMAKE_BIN="${SNAKEMAKE_BIN:-snakemake}"
SNAKEMAKE_CORES="${SNAKEMAKE_CORES:-1}"
exec "$SNAKEMAKE_BIN" \
  --snakefile "$PROJECT_ROOT/workflow/Snakefile" \
  --directory "$PROJECT_ROOT" \
  --configfile "$PROJECT_ROOT/config/config.yaml" \
  --cores "$SNAKEMAKE_CORES" \
  --until collectri \
  "$@"
