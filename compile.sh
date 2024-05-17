#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: $0 [Debug|Release]"
    exit 1
}

# Check if build directory exists
if [ -d "build" ]; then
    build_exists=true
else
    build_exists=false
fi

# Check if build directory exists and delete if it does (only if an argument is provided)
if [ -n "$1" ] && [ -d "build" ]; then
    echo "Removing existing build directory..."
    rm -rf build || { echo "Failed to remove build directory"; exit 1; }
    build_exists=false
fi

# Create build directory if it doesn't exist and if a build folder didn't exist before
if ! $build_exists; then
    echo "Create build directory..."
    mkdir -p build || { echo "Failed to create build directory"; exit 1; }
    build_exists=true
fi

# Determine build type
if [ -n "$1" ]; then
    case "$1" in
        Debug)
            build_type="Debug"
            ;;
        Release)
            build_type="Release"
            ;;
        *)
            echo "Invalid build type. Please specify Debug or Release."
            usage
            ;;
    esac
else
    if [ "$build_exists" != true ]; then
        echo "No build type specified. Please provide either Debug or Release."
        usage
    fi
fi

if [ -n "$build_type" ]; then
    echo "Configuring with $build_type build type..."
    cmake -S ${PWD} -B ${PWD}/build -DCMAKE_BUILD_TYPE="$build_type" || { echo "CMake configuration failed"; exit 1; }
fi

cmake --build ${PWD}/build -j 16



# Configure cmake
#cmake -DCMAKE_BUILD_TYPE=$build_type ..

# # Check if cmake configuration was successful
# if [ $? -ne 0 ]; then
#     echo "Failed to configure cmake"
#     exit 1
# fi

# # Build the project
# cmake --build . --config $build_type

# # Check if build was successful
# if [ $? -ne 0 ]; then
#     echo "Build failed"
#     exit 1
# fi

# # Optionally, run tests or additional commands here

# echo "Build completed successfully."



# # Configure with specified build type if provided
# if [ -n "$build_type" ]; then
#     echo "Configuring with $build_type build type..."
#     cmake -S . -B build -DCMAKE_BUILD_TYPE="$build_type" || { echo "CMake configuration failed"; exit 1; }
# fi

# # Build the project
# echo "Building the project..."
# # Use the default number of parallel jobs
# cmake --build build --parallel || { echo "Build failed"; exit 1; }




# #!/bin/bash

# # Function to display usage information
# usage() {
#     echo "Usage: $0 [Debug|Release]"
#     exit 1
# }

# # Check if build directory exists
# if [ -d "build" ]; then
#     build_exists=true
# else
#     build_exists=false
# fi

# # Check if build directory exists and delete if it does (only if an argument is provided)
# if [ -n "$1" ] && [ -d "build" ]; then
#     echo "Removing existing build directory..."
#     rm -rf build || { echo "Failed to remove build directory"; exit 1; }
# fi

# # Create build directory if it doesn't exist and if a build folder didn't exist before
# if ! $build_exists; then
#     mkdir -p build || { echo "Failed to create build directory"; exit 1; }
# fi

# # Determine build type
# if [ -n "$1" ]; then
#     case "$1" in
#         Debug)
#             build_type="Debug"
#             ;;
#         Release)
#             build_type="Release"
#             ;;
#         *)
#             echo "Invalid build type. Please specify Debug or Release."
#             usage
#             ;;
#     esac
# else
#     if [ "$build_exists" != true ]; then
#         echo "No build type specified. Please provide either Debug or Release."
#         usage
#     fi
# fi

# # Configure with specified build type if provided
# # if [ -n "$build_type" ]; then
# #     echo "Configuring with $build_type build type..."
# #     cmake -S . -B build -DCMAKE_BUILD_TYPE="$build_type" -- -j16 || { echo "CMake configuration failed"; exit 1; }
# # fi

# # Configure with specified build type if provided
# if [ -n "$build_type" ]; then
#     echo "Configuring with $build_type build type..."
#     cmake -S . -B build -DCMAKE_BUILD_TYPE="$build_type" || { echo "CMake configuration failed"; exit 1; }
# fi

# # Build the projec
# echo "Building the project..."
# cmake --build build || { echo "Build failed"; exit 1; }
