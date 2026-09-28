#!/usr/bin/env bash
# Portable wrapper; each run writes into a new staging directory.
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONIOENCODING=utf-8
exec python "$SCRIPT_DIR/regenerate_chi4.py" "$@"
