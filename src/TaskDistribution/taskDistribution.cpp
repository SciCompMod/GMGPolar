#include "../../include/TaskDistribution/taskDistribution.h"

// Generates a uniform Task Distribution for OpenMP multithreading.

TaskDistribution::TaskDistribution(){}

TaskDistribution::TaskDistribution(const int numTasks, const int minChunkSize, const int numThreads)
{
    partition_.resize(numThreads + 1);

    int chunkSize = numTasks / numThreads;
    int current = 0;
    partition_[0] = 0;
    for (int i = 0; i < numThreads; ++i) {
        current += std::max(chunkSize, minChunkSize);
        partition_[i+1] = current;
    }

    int threadsUsed = numThreads;
    for (int i = 0; i < numThreads+1; i ++) {
        if (partition_[i] <= numTasks)
            threadsUsed = i;
        else
            break;
    }

    int increasing = numTasks - partition_[threadsUsed];

    for (int i = 0; i < increasing; i++){
        for (int j = std::min(i+1, threadsUsed); j < threadsUsed+1; j++){
            partition_[j] += 1;
        }
    }

    for (int i = threadsUsed+1; i < numThreads+1; i++){
        partition_[i] = numTasks;
    }

    partition_[0] = 0;
}

int TaskDistribution::getStart(const int threadID) const {
    assert(threadID >= 0 && static_cast<size_t>(threadID) < partition_.size()-1);

    return partition_[threadID];
}

int TaskDistribution::getEnd(const int threadID) const {
    assert(threadID >= 0 && static_cast<size_t>(threadID) < partition_.size()-1);

    return partition_[threadID+1];
}

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