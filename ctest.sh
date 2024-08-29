#!/bin/bash

# Here we test for rare race conditions!

# Change to the build directory
cd build

# Loop to run ctest 100 times
for i in {1..100}
do
    echo "Running test iteration $i"
    ctest

    # Check if ctest failed
    if [ $? -ne 0 ]; then
        echo "ctest failed on iteration $i"
        cd ..
        exit 1
    fi
done

# Return to the previous directory
cd ..

echo "All 100 tests passed successfully."