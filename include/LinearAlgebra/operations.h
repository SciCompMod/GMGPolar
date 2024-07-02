#pragma once

#include <cassert>
#include <omp.h>

#include "vector.h"
#include "matrix.h"
#include "../common/equals.h"

#include "../TaskDistribution/taskDistribution.h"

template<typename T>
bool equals(const Vector<T>& lhs, const Vector<T>& rhs){
    const int size_l = lhs.size();
    const int size_r = rhs.size();
    if(size_l != size_r){
        return false;
    }
    bool result = true;
    #pragma omp parallel for shared(result)
    for(int i = 0; i < size_l; ++i) {
        if (!equals(lhs[i], rhs[i])) {
            #pragma omp atomic write
            result = false;
        }
    }
    return result;
}

template<typename T>
void assign(Vector<T>& lhs, const T& value) {
    const int maxThreads = omp_get_max_threads();
    const int minimalChunkSize = 64;
    
    TaskDistribution partition(lhs.size(), minimalChunkSize, maxThreads);
    
    #pragma omp parallel num_threads(maxThreads)
    {
        const int threadID = omp_get_thread_num();
        const int local_start = partition.getStart(threadID);
        const int local_end = partition.getEnd(threadID);
        
        for (int i = local_start; i < local_end; ++i) {
            lhs[i] = value;
        }
    }
}

template<typename T>
void add(Vector<T>& result, const Vector<T>& x){
    assert(result.size() == x.size());

    const int maxThreads = omp_get_max_threads();
    const int minimalChunkSize = 64;
    
    TaskDistribution partition(result.size(), minimalChunkSize, maxThreads);
    
    #pragma omp parallel num_threads(maxThreads)
    {
        const int threadID = omp_get_thread_num();
        const int local_start = partition.getStart(threadID);
        const int local_end = partition.getEnd(threadID);
        
        for (int i = local_start; i < local_end; ++i) {
            result[i] += x[i];
        }
    }
}

template<typename T>
void subtract(Vector<T>& result, const Vector<T>& x){
    assert(result.size() == x.size());

    const int maxThreads = omp_get_max_threads();
    const int minimalChunkSize = 64;
    
    TaskDistribution partition(result.size(), minimalChunkSize, maxThreads);
    
    #pragma omp parallel num_threads(maxThreads)
    {
        const int threadID = omp_get_thread_num();
        const int local_start = partition.getStart(threadID);
        const int local_end = partition.getEnd(threadID);
        
        for (int i = local_start; i < local_end; ++i) {
            result[i] -= x[i];
        }
    }
}


/* x = alpha * x + beta * y */
template<typename T>
void linear_combination(Vector<T>& x, const T& alpha, const Vector<T>& y, const T& beta){
    assert(x.size() == y.size());

    const int maxThreads = omp_get_max_threads();
    const int minimalChunkSize = 64;
    
    TaskDistribution partition(x.size(), minimalChunkSize, maxThreads);
    
    #pragma omp parallel num_threads(maxThreads)
    {
        const int threadID = omp_get_thread_num();
        const int local_start = partition.getStart(threadID);
        const int local_end = partition.getEnd(threadID);
        
        for (int i = local_start; i < local_end; ++i) {
            x[i] *= alpha;
            x[i] += beta * y[i];
        }
    }
}

/* x = alpha * x */
template<typename T>
void multiply(Vector<T>& x, const T& alpha){
    const int maxThreads = omp_get_max_threads();
    const int minimalChunkSize = 64;
    
    TaskDistribution partition(x.size(), minimalChunkSize, maxThreads);
    
    #pragma omp parallel num_threads(maxThreads)
    {
        const int threadID = omp_get_thread_num();
        const int local_start = partition.getStart(threadID);
        const int local_end = partition.getEnd(threadID);
        
        for (int i = local_start; i < local_end; ++i) {
            x[i] *= alpha;
        }
    }
}

template<typename T>
T dot_product(const Vector<T>& lhs, const Vector<T>& rhs){
    assert(lhs.size() == rhs.size());

    const int maxThreads = omp_get_max_threads();
    const int minimalChunkSize = 64;
    
    TaskDistribution partition(lhs.size(), minimalChunkSize, maxThreads);

    T global_sum = 0.0;
    
    #pragma omp parallel num_threads(maxThreads)
    {
        const int threadID = omp_get_thread_num();
        const int local_start = partition.getStart(threadID);
        const int local_end = partition.getEnd(threadID);

        T local_sum = 0.0;
        
        for (int i = local_start; i < local_end; ++i) {
            local_sum += lhs[i] * rhs[i];
        }

        #pragma omp critical
        {
            global_sum += local_sum;
        }
    }
    return global_sum;
}

template<typename T>
T l1_norm(const Vector<T>& x){
    const int maxThreads = omp_get_max_threads();
    const int minimalChunkSize = 64;
    
    TaskDistribution partition(x.size(), minimalChunkSize, maxThreads);

    T global_sum = 0.0;
    
    #pragma omp parallel num_threads(maxThreads)
    {
        const int threadID = omp_get_thread_num();
        const int local_start = partition.getStart(threadID);
        const int local_end = partition.getEnd(threadID);

        T local_sum = 0.0;
        
        for (int i = local_start; i < local_end; ++i) {
            local_sum += fabs(x[i]);
        }

        #pragma omp critical
        {
            global_sum += local_sum;
        }
    }
    return global_sum;
}

template<typename T>
T l2_norm_squared(const Vector<T>& x){
    const int maxThreads = omp_get_max_threads();
    const int minimalChunkSize = 64;
    
    TaskDistribution partition(x.size(), minimalChunkSize, maxThreads);

    T global_sum = 0.0;
    
    #pragma omp parallel num_threads(maxThreads)
    {
        const int threadID = omp_get_thread_num();
        const int local_start = partition.getStart(threadID);
        const int local_end = partition.getEnd(threadID);

        T local_sum = 0.0;
        
        for (int i = local_start; i < local_end; ++i) {
            local_sum += x[i] * x[i];
        }

        #pragma omp critical
        {
            global_sum += local_sum;
        }
    }
    return global_sum;
}

template<typename T>
T infinity_norm(const Vector<T>& x){
    const int maxThreads = omp_get_max_threads();
    const int minimalChunkSize = 64;
    
    TaskDistribution partition(x.size(), minimalChunkSize, maxThreads);

    T global_max = 0.0;
    
    #pragma omp parallel num_threads(maxThreads)
    {
        const int threadID = omp_get_thread_num();
        const int local_start = partition.getStart(threadID);
        const int local_end = partition.getEnd(threadID);

        T local_max = 0.0;
        
        for (int i = local_start; i < local_end; ++i) {
            T absolute_value = fabs(x[i]);
            if(local_max < absolute_value) local_max = absolute_value;
        }

        #pragma omp critical
        {
            if(global_max < local_max) global_max = local_max;
        }
    }
    return global_max;
}



template<typename T>
void multiply(Vector<T>& result, const SparseMatrix<T>& lhs, const Vector<T>& rhs){
    for (int i = 0; i < lhs.non_zero_size(); i++){
        result[lhs.row_index(i)-1] += lhs.value(i) * rhs[lhs.col_index(i)-1];
    }
}