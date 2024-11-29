#!/bin/bash

# Navigate to the build directory
cd ../build || { echo "Build folder not found!"; exit 1; }

# Run ctest in verbose mode
ctest --verbose