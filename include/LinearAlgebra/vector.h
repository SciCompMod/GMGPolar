#pragma once

#include <algorithm>
#include <cassert>
#include <initializer_list>
#include <iostream>
#include <memory>
#include <vector>

template <typename T>
class Vector
{
public:
    Vector();
    Vector(const Vector& other);
    Vector(Vector&& other) noexcept;

    explicit Vector(std::size_t size);
    Vector(const std::initializer_list<T>& init);
    Vector(const std::vector<T>& init);

    Vector& operator=(const Vector& other);
    Vector& operator=(Vector&& other) noexcept;

    ~Vector() = default;

    std::size_t size() const noexcept;

    inline const T& operator[](std::size_t index) const;
    inline T& operator[](std::size_t index);

    T* begin() noexcept;
    T* end() noexcept;
    const T* begin() const noexcept;
    const T* end() const noexcept;

    template <typename U>
    friend std::ostream& operator<<(std::ostream& stream, const Vector<U>& vector);

private:
    std::size_t size_;
    std::unique_ptr<T[]> values_;
};

template <typename U>
std::ostream& operator<<(std::ostream& stream, const Vector<U>& vector)
{
    stream << "[";
    for (std::size_t i = 0; i < vector.size(); ++i) {
        stream << vector[i];
        if (i + 1 != vector.size())
            stream << ", ";
    }
    stream << "]";
    return stream;
}

// default construction
template <typename T>
Vector<T>::Vector()
    : size_(0)
    , values_(std::make_unique<T[]>(0))
{
}

// copy construction
template <typename T>
Vector<T>::Vector(const Vector& other)
    : size_(other.size_)
    , values_(std::make_unique<T[]>(size_))
{
#pragma omp parallel for if (size_ > 10'000)
    for (std::size_t i = 0; i < size_; ++i) {
        values_[i] = other.values_[i];
    }
}

// copy assignment
template <typename T>
Vector<T>& Vector<T>::operator=(const Vector& other)
{
    if (this == &other) {
        return *this;
    }

    if (size_ != other.size_) {
        size_   = other.size_;
        values_ = std::make_unique<T[]>(size_);
    }

#pragma omp parallel for if (size_ > 10'000)
    for (std::size_t i = 0; i < size_; ++i) {
        values_[i] = other.values_[i];
    }

    return *this;
}

// move construction
template <typename T>
Vector<T>::Vector(Vector&& other) noexcept
    : size_(other.size_)
    , values_(std::move(other.values_))
{
    other.size_ = 0;
}

// move assignment
template <typename T>
Vector<T>& Vector<T>::operator=(Vector&& other) noexcept
{
    if (this != &other) {
        size_       = other.size_;
        values_     = std::move(other.values_);
        other.size_ = 0;
    }
    return *this;
}

// constructor with size
template <typename T>
Vector<T>::Vector(std::size_t size)
    : size_(size)
    , values_(std::make_unique<T[]>(size_))
{
}

// constructor with initializer list
template <typename T>
Vector<T>::Vector(const std::initializer_list<T>& init)
    : size_(init.size())
    , values_(std::make_unique<T[]>(size_))
{
    std::copy(init.begin(), init.end(), values_.get());
}

// constructor with std::vector
template <typename T>
Vector<T>::Vector(const std::vector<T>& init)
    : size_(init.size())
    , values_(std::make_unique<T[]>(size_))
{
    std::copy(init.begin(), init.end(), values_.get());
}

// getter []
template <typename T>
inline const T& Vector<T>::operator[](std::size_t index) const
{
    assert(index < size_);
    return values_[index];
}

// setter []
template <typename T>
inline T& Vector<T>::operator[](std::size_t index)
{
    assert(index < size_);
    return values_[index];
}

// get vector's size
template <typename T>
std::size_t Vector<T>::size() const noexcept
{
    return size_;
}

template <typename T>
T* Vector<T>::begin() noexcept
{
    return values_.get();
}
template <typename T>
T* Vector<T>::end() noexcept
{
    return values_.get() + size_;
}

template <typename T>
const T* Vector<T>::begin() const noexcept
{
    return values_.get();
}
template <typename T>
const T* Vector<T>::end() const noexcept
{
    return values_.get() + size_;
}
