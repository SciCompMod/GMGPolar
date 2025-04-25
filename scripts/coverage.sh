#!/bin/bash

# Set up environment variables for dependencies
# export MUMPS_DIR="$HOME/spack/opt/spack/linux-ubuntu24.04-x86_64/gcc-13.3.0/mumps-5.5.1-2v3iw45om32klbzpdsssblfbm2zklrl6"
# export METIS_DIR="$HOME/spack/opt/spack/linux-ubuntu24.04-x86_64_v4/gcc-13.3.0/metis-5.1.0-ab7tqwpnvjmngmrgxodnjykdsysypfbx"
# export LIKWID_DIR="$HOME/spack/opt/spack/linux-ubuntu24.04-x86_64_v4/gcc-13.3.0/likwid-5.2.2-nswn2ksoi2phc7qx55urkko6hxvfqkhf"

# Default values
CLEAN=True
BUILD_TYPE="Debug"
JOBS=16
TESTS=True
COVERAGE=True
MUMPS=False
LIKWID=False

# Validate boolean options
for var in CLEAN COVERAGE LIKWID MUMPS TESTS; do
    declare -n value=$var
    if [[ "$value" != "True" && "$value" != "False" ]]; then
        echo "Invalid value for $var: $value (must be True or False)"
        exit 1
    fi
done

# Validate build type
case "$BUILD_TYPE" in
    Debug|Release)
        ;;
    *)
        echo "Invalid build type: $BUILD_TYPE"
        usage
        ;;
esac

# Set up build directory
BUILD_DIR="../build"

if [ "$CLEAN" = "True" ] && [ -d "$BUILD_DIR" ]; then
    echo "Removing existing build directory..."
    rm -rf "$BUILD_DIR" || { echo "Failed to remove build directory"; exit 1; }
fi

if [ ! -d "$BUILD_DIR" ]; then
    echo "Creating build directory..."
    mkdir -p "$BUILD_DIR" || { echo "Failed to create build directory"; exit 1; }
fi

# Configure with CMake
echo "Configuring with options:"
echo "  Build type: $BUILD_TYPE"
echo "  Tests: $TESTS"
echo "  Coverage: $COVERAGE"
echo "  LIKWID: $LIKWID"
echo "  MUMPS: $MUMPS"

CMAKE_OPTIONS=(
    "-DCMAKE_BUILD_TYPE=$BUILD_TYPE"
    "-DGMGPOLAR_BUILD_TESTS=$TESTS"
    "-DGMGPOLAR_USE_LIKWID=$LIKWID"
    "-DGMGPOLAR_USE_MUMPS=$MUMPS"
    "-DGMGPOLAR_ENABLE_COVERAGE=$COVERAGE"
)

cmake -S .. -B "$BUILD_DIR" "${CMAKE_OPTIONS[@]}" || { echo "CMake configuration failed"; exit 1; }

# Build the project
echo "Building with $JOBS parallel jobs..."
cmake --build "$BUILD_DIR" -j "$JOBS" || { echo "Build failed"; exit 1; }

# Run coverage if enabled
if [ "$COVERAGE" = "True" ] && [ "$TESTS" = "True" ]; then
    echo "Generating coverage report..."
    cmake --build "$BUILD_DIR" --target coverage || { echo "Coverage generation failed"; exit 1; }
    
    # Open the coverage report if possible
    if command -v xdg-open >/dev/null 2>&1; then
        xdg-open "$BUILD_DIR/coverage-report/index.html"
    else
        echo "Coverage report generated at: $BUILD_DIR/coverage-report/index.html"
    fi
elif [ "$COVERAGE" = "True" ] && [ "$TESTS" = "False" ]; then
    echo "Warning: Coverage requested but tests are disabled"
fi

echo "Build completed successfully!"