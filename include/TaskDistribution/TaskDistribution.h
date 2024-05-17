#pragma once

#include <omp.h>
#include <vector>
#include <cassert>
#include <cstddef>
#include <iostream>

class TaskDistribution {
public:
    TaskDistribution();
    TaskDistribution(const int numTasks, const int minChunkSize, const int numThreads);
    int getStart(const int threadID) const;
    int getEnd(const int threadID) const;

private:
    std::vector<int> partition_;
};

// Generates a uniform Task Distribution for OpenMP multithreading.

// Example 1. With numTasks = 10, minChunkSize = 1 and numThreads = 4 we get
// Thread 0: start = 0. end = 3
// Thread 1: start = 3. end = 6
// Thread 2: start = 6. end = 8
// Thread 3: start = 8. end = 10
// The size of the partitions is then [3,3,2,2].

// Example 2. With numTasks = 10, minChunkSize = 3 and numThreads = 4 we get
// Thread 0: start = 0. end = 4
// Thread 1: start = 4. end = 7
// Thread 2: start = 7. end = 10
// Thread 3: start = 10. end = 10
// The size of the partitions is then [4,3,3,0].

// Example 3. With numTasks = 10, minChunkSize = 4 and numThreads = 4 we get
// Thread 0: start = 0. end = 5
// Thread 1: start = 5. end = 10
// Thread 2: start = 10. end = 10
// Thread 3: start = 10. end = 10
// The size of the partitions is then [5,5,0,0].

// Example 4. With numTasks = 10, minChunkSize = 15 and numThreads = 4 we get
// Thread 0: start = 0. end = 10
// Thread 1: start = 10. end = 10
// Thread 2: start = 10. end = 10
// Thread 3: start = 10. end = 10
// The size of the partitions is then [10,0,0,0].

// Consider using for(int i = start, i < end; i++){ ... }

// Note: If numTasks is less than minChunkSize,
// the first thread does not fulfill the minimal chunk size requirement.