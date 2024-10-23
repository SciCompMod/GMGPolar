#include "../../include/PolarGrid/point.h"

Point::Point(double value)
{
    assert(space_dimension >= 0);
    std::fill(std::begin(data_), std::end(data_), value);
}

Point::Point(double x, double y)
{
    assert(space_dimension == 2);
    data_[0] = x;
    data_[1] = y;
}

Point::Point(double x, double y, double z)
{
    assert(space_dimension == 3);
    data_[0] = x;
    data_[1] = y;
    data_[2] = z;
}

int Point::size() const
{
    return space_dimension;
}

const double& Point::operator[](int i) const
{
    assert(i >= 0);
    assert(i < space_dimension);
    return data_[i];
}

double& Point::operator[](int i)
{
    assert(i >= 0);
    assert(i < space_dimension);
    return data_[i];
}

bool equals(const Point& lhs, const Point& rhs)
{
    for (int i = 0; i < space_dimension; ++i)
    {
        if (!equals(lhs[i], rhs[i])) return false;
    }
    return true;
}

double norm(const Point& point)
{
    double result = 0;
    for (int i = 0; i < space_dimension; ++i)
    {
        result += point[i] * point[i];
    }
    return sqrt(result);
}

void add(Point& result, const Point& lhs, const Point& rhs)
{
    assert(lhs.size() == rhs.size());
    assert(result.size() == lhs.size());
    for (int i = 0; i < space_dimension; i++)
    {
        result[i] = lhs[i] + rhs[i];
    }
}

void subtract(Point& result, const Point& lhs, const Point& rhs)
{
    assert(lhs.size() == rhs.size());
    assert(result.size() == lhs.size());
    for (int i = 0; i < space_dimension; i++)
    {
        result[i] = lhs[i] - rhs[i];
    }
}

void multiply(Point& result, const double& scalar, const Point& lhs)
{
    assert(result.size() == lhs.size());
    for (int i = 0; i < space_dimension; i++)
    {
        result[i] = scalar * lhs[i];
    }
}