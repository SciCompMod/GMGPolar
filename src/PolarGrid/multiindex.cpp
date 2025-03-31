#include "../../include/PolarGrid/multiindex.h"

MultiIndex::MultiIndex(int value)
{
    assert(space_dimension >= 0);
    std::fill(std::begin(data_), std::end(data_), value);
}

MultiIndex::MultiIndex(int i, int j)
{
    assert(space_dimension == 2);
    data_[0] = i;
    data_[1] = j;
}

MultiIndex::MultiIndex(int i, int j, int k)
{
    assert(space_dimension == 3);
    data_[0] = i;
    data_[1] = j;
    data_[2] = k;
}

int MultiIndex::size() const
{
    return space_dimension;
}

bool MultiIndex::operator==(const MultiIndex& other) const
{
    for (int d = 0; d < space_dimension; d++) {
        if (data_[d] != other.data_[d]) {
            return false;
        }
    }
    return true;
}

bool MultiIndex::operator!=(const MultiIndex& other) const
{
    for (int d = 0; d < space_dimension; d++) {
        if (data_[d] != other.data_[d]) {
            return true;
        }
    }
    return false;
}

const int& MultiIndex::operator[](int i) const
{
    assert(i >= 0);
    assert(i < space_dimension);
    return data_[i];
}

int& MultiIndex::operator[](int i)
{
    assert(i >= 0);
    assert(i < space_dimension);
    return data_[i];
}

int* MultiIndex::data()
{
    return data_;
}

const int* MultiIndex::data() const
{
    return data_;
}