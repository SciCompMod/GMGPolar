#include "../../include/PolarGrid/polargrid.h"

void PolarGrid::writeToFile(const std::string& file_r, const std::string& file_theta, const int precision) const
{
    writeVectorToFile(file_r, radii_, precision);
    writeVectorToFile(file_theta, angles_, precision);
}

void PolarGrid::writeVectorToFile(const std::string& filename, const std::vector<double>& vector, const int precision) const
{
    // Open the file for writing
    std::ofstream outputFile(filename);

    // Check if the file is opened successfully
    if (!outputFile.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Set the precision for output file
    outputFile << std::fixed << std::setprecision(precision);

    // Write each double from the vector to the file
    for (const auto& num : vector)
    {
        outputFile << num << std::endl;
    }

    // Close the file
    outputFile.close();
}

void PolarGrid::loadVectorFromFile(const std::string& filename, std::vector<double>& vector) const
{
    // Open the file for reading
    std::ifstream inputFile(filename);

    // Check if the file is opened successfully
    if (!inputFile.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Clear the vector before loading new data
    vector.clear();

    // Read data from the file into the vector
    double value;
    while (inputFile >> value)
    {
        vector.push_back(value);
    }

    // Close the file
    inputFile.close();
}