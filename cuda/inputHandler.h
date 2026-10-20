#ifndef INPUTHANDLER_H
#define INPUTHANDLER_H

#include <string>

class InputValues {

public:
    int numSideNodes;
    int maxItt;
    double initialDensity;
    double reynoldsNumber;
    double tau;
    double maxError;

    InputValues(const std::string& inputFileName);
private:
    void readInputFile(const std::string& inputFileName);

};

#endif
