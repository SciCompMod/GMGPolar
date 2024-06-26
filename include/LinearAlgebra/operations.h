#pragma once

#include <cassert>
#include <omp.h>

#include "vector.h"
#include "matrix.h"
#include "../common/equals.h"

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
void assign(Vector<T>& lhs, const T& rhs){
    const int size = lhs.size();
    #pragma omp parallel for
    for (int i = 0; i < size; i++){
        lhs[i]=rhs;
    }
}

template<typename T>
void add(Vector<T>& result, const Vector<T>& x){
    assert(result.size()==x.size());
    const int size=x.size();
    #pragma omp parallel for
    for(int i = 0; i < size; ++i) {
        result[i] += x[i];
    }
}

template<typename T>
void add(Vector<T>& result, const Vector<T>& lhs,const Vector<T>& rhs){
    assert(lhs.size()==rhs.size());
    assert(result.size()==rhs.size());
    const int size=lhs.size();
    #pragma omp parallel for
    for(int i = 0; i < size; ++i) {
        result[i] = lhs[i] + rhs[i];
    }
}

template<typename T>
void subtract(Vector<T>& result, const Vector<T>& x){
    assert(result.size()==x.size());
    const int size=x.size();
    #pragma omp parallel for
    for(int i=0;i<size;++i){
        result[i] -= x[i];
    }
}

template<typename T>
void subtract(Vector<T>& result, const Vector<T>& lhs,const Vector<T>& rhs){
    assert(lhs.size()==rhs.size());
    assert(result.size()==rhs.size());
    const int size=lhs.size();
    #pragma omp parallel for
    for(int i=0;i<size;++i){
        result[i] = lhs[i] - rhs[i];
    }
}

template<typename T>
T dot_product(const Vector<T>& lhs, const Vector<T>& rhs){
    assert(lhs.size() == rhs.size());
    T sum = 0;
    const int size = lhs.size();
    #pragma omp parallel for reduction(+:sum)
    for(int i = 0; i < size; ++i) {
        sum += lhs[i] * rhs[i];
    }
    return sum;
}

template<typename T>
void multiply(Vector<T>& result, const Vector<T>& lhs, const T& rhs){
    assert(result.size()==lhs.size());
    const int size=lhs.size();
    #pragma omp parallel for
    for(int i=0;i<size;++i){
        result[i]=rhs*lhs[i];
    }
}


template<typename T>
void multiply(Vector<T>& result, const SparseMatrix<T>& lhs, const Vector<T>& rhs){
    assert(result.size()==lhs.rows());
    assert(lhs.size()()==rhs.cols());
    for (int i = 0; i < lhs.non_zero_size(); i++)
    {
        result[lhs.row_index(i)-1] += lhs.value(i) * rhs[lhs.col_index(i)-1];
    }
}

