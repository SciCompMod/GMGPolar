#pragma once

#include <memory>
#include <cassert>
#include <algorithm>
#include <initializer_list>
#include <sstream>
#include <vector>
#include <iostream>

template<typename T>
class Vector
{
public:
    Vector();
    Vector(const Vector& other);
    Vector(Vector&& other) noexcept;

    explicit Vector(int size);
    Vector(const std::initializer_list<T>& init);
    Vector(const std::vector<T>& init);

    Vector& operator=(const Vector& other);
    Vector& operator=(Vector&& other) noexcept;

    int size() const noexcept;

    inline const T& operator[](int index) const;
    inline T& operator[](int index);

    T* begin() noexcept;
    T* end() noexcept;
    const T* begin() const noexcept;
    const T* end() const noexcept;

    template<typename U>
    friend std::ostream& operator<<(std::ostream& stream, const Vector<U>& vector);

private:
    int size_;
    std::unique_ptr<T[]> values_;
};

template<typename U>
std::ostream& operator<<(std::ostream& stream, const Vector<U>& vector)
{
    stream << "[";
    for (int i = 0; i < vector.size(); ++i) { 
        stream << vector[i]; 
        if (i != vector.size() - 1) 
            stream << ", "; 
    }
    stream << "]\n";
    return stream;
}

// default construction
template<typename T>
Vector<T>::Vector():
    size_(0),
    values_(nullptr)
{}

// copy construction
template<typename T>
Vector<T>::Vector(const Vector& other):
    size_(other.size_),
    values_(std::make_unique<T[]>(size_))
{
    #pragma omp parallel for if(size_ > 100'000)
    for (std::size_t i = 0; i < size_; ++i) {
        values_[i] = other.values_[i];
    }
}

// copy assignment
template<typename T>
Vector<T>& Vector<T>::operator=(const Vector& other) {
    if (this == &other) {
        /* Self-assignment, no work needed */
        return *this;
    }
    
    /* Allocate new memory if necessary */
    if (size_ != other.size_) {
        size_ = other.size_;
        values_ = std::make_unique<T[]>(size_);
    }
    
    #pragma omp parallel for if(size_ > 100'000)
    for (std::size_t i = 0; i < size_; ++i) {
        values_[i] = other.values_[i];
    }

    return *this;
}

// move construction
template<typename T>
Vector<T>::Vector(Vector&& other) noexcept:
    size_(other.size_),
    values_(std::move(other.values_))
{
    other.size_ = 0;
}

// move assignment
template<typename T>
Vector<T>& Vector<T>::operator=(Vector&& other) noexcept {
    if (this != &other) {
        size_ = other.size_;
        values_ = std::move(other.values_);
        other.size_ = 0;
    }
    return *this;
}

// constrcutor with size
template<typename T>
Vector<T>::Vector(int size):
    size_(size),
    values_(std::make_unique<T[]>(size_))
{}

// constructor with initialization
template<typename T>
Vector<T>::Vector(const std::initializer_list<T>& init):
    size_(static_cast<int>(init.size())),
    values_(std::make_unique<T[]>(size_))
{
    std::copy(init.begin(), init.end(), values_.get());
}
// constructor with std::vector initialization
template<typename T>
Vector<T>::Vector(const std::vector<T>& init):
    size_(static_cast<int>(init.size())),
    values_(std::make_unique<T[]>(size_))
{
    assert(!init.empty());
    std::copy(init.begin(), init.end(), values_.get());
}

// getter []
template<typename T>
inline const T& Vector<T>::operator[](int index) const{
    assert(index >= 0); assert(index < size_);
    return values_[index];
}
// setter []
template<typename T>
inline T& Vector<T>::operator[](int index){
    assert(index >= 0); assert(index < size_);
    return values_[index];
}
// get vector's size
template<typename T>
int Vector<T>::size() const noexcept{ 
    return size_;
}

template<typename T>
T* Vector<T>::begin() noexcept {
    return values_.get();
}
template<typename T>
T* Vector<T>::end() noexcept {
    return values_.get() + size_;
}

template<typename T>
const  T* Vector<T>::begin() const noexcept {
    return values_.get();
}

template<typename T>
const T* Vector<T>::end() const noexcept {
    return values_.get() + size_;
}