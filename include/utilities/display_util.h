#pragma once

#include <iostream>
#include <string>
#include <vector>

/***************************************************************************
* Display arrays or vectors
**************************************************************************/
/*!
 *  \brief Forwards all values inside an array to std::cout.
 * 
 * \param na: Size of the array.
 * \param a: The array.
 * \param s_a: Name of the array (to be printed).
 *
 */
template <typename T>
void disp(int na, T* a, const std::string& s_a)
{
    std::cout << s_a << "(" << na << "): ";
    for (int i = 0; i < na; i++)
        std::cout << a[i] << " ";
    std::cout << "\n";
} 

/*!
 *  \brief Forwards all values inside a vector to std::cout.
 * 
 * \param a: The vector.
 * \param s_a: Name of the vector (to be printed).
 *
 */
template <typename T>
void disp(std::vector<T> a, const std::string& s_a)
{
    std::cout << s_a << "(" << a.size() << "): ";
    for (std::size_t i = 0; i < a.size(); i++)
        std::cout << a[i] << " ";
    std::cout << "\n";
} 
