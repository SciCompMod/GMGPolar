#!/bin/bash

# Get the directory of the script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Set path to the GMGPolar executable
BUILD_DIRECTORY="${SCRIPT_DIR}/../../build"

# Adjust parameters in src/convergence_order.cpp

"$BUILD_DIRECTORY"/convergence_order