#!/bin/bash

# Adjust parameters in src/convergence_order.cpp

# Get the directory of the script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Set path to the GMGPolar executable
CONVERGENCE_EXEC="${SCRIPT_DIR}/../../build/convergence_order"

"$CONVERGENCE_EXEC"