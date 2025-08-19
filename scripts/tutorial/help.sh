#!/bin/bash

# Get the directory of the script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Set path to the GMGPolar executable
HELP_EXEC="${SCRIPT_DIR}/../../build/gmgpolar"

"$HELP_EXEC" --help