#pragma once

#include <algorithm>
#include <cassert>
#include <iterator>

#include "../common/equals.h"
#include "../common/space_dimension.h"

class Point
{
public:
    static_assert(space_dimension > 0 && space_dimension <= 3, "Invalid space dimension");

    //! Creates a point with undefined coordinate values.
    Point() = default;

    //! Initializes all coordinate values of the point with `value`.
    explicit Point(double value);
    //! Initializes 2-dimensional point coordinate values with `x` and `y`.
    //! Works only if `space_dimension` is 2.
    explicit Point(double x, double y);
    //! Initializes 3-dimensional point coordinate values with `x`, `y` and `z`.
    //! Works only if `space_dimension` is 3.
    explicit Point(double x, double y, double z);

    //! Returns the size of the point.
    //! This is equal to `space_dimension`.
    int size() const;

    //! Returns the `i`th coordinate value of the point.
    double operator[](int i) const;

    //! Returns a mutable reference to the `i`th coordinate value of the point.
    double& operator[](int i);

private:
    double data_[space_dimension];
};

bool equals(const Point& lhs, const Point& rhs);
double norm(const Point& point);
void add(Point& result, const Point& lhs, const Point& rhs);
void subtract(Point& result, const Point& lhs, const Point& rhs);
void multiply(Point& result, double scalar, const Point& lhs);