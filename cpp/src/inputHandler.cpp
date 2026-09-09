#include "inputHandler.h"

#include <fstream>
#include <vector>
#include <stdexcept>

InputValues::InputValues(const std::string& inputFileName) {
    readInputFile(inputFileName);
}

void InputValues::readInputFile(const std::string& inputFileName) {
    std::string inputFilePath = "../inputs/"+inputFileName+".txt";

    std::ifstream fileReader(inputFilePath);

    if (!fileReader.is_open()) {
        throw std::runtime_error("Could not open input file: " + inputFilePath);
    }

    std::vector<std::string> fileData;

    std::string line = "";
    std::string lineData = "";
    std::size_t leftDelimiter;
    std::size_t rightDelimiter;
    while (std::getline(fileReader, line)) {
        leftDelimiter = line.find('[');
        rightDelimiter = line.find(']');
        lineData = line.substr(leftDelimiter+1, rightDelimiter-leftDelimiter-1);

        fileData.push_back(lineData);
    }

    fileReader.close();

    numSideNodes = stoi(fileData.at(0));
    maxItt = stoi(fileData.at(1));
    initialDensity = stod(fileData.at(2));
    reynoldsNumber = stod(fileData.at(3));
    tau = stod(fileData.at(4));
    maxError = stod(fileData.at(5));

    return;
}
