#include "inputHandler.h"
#include "matrixMethods/matrix.h"

#include <string>
#include <iostream>

int main () {

    std::cout << "MAIN STARTED -----\n";

    std::string inputFileName = "case0";

    InputValues inputs = InputValues(inputFileName);

    // Read Inputs
    // note to self: change to proper notation if const doesnt cause issues
    const int numX = inputs.numSideNodes;
    const int numY = numX;
    const int maxItt = inputs.maxItt;
    const int initialDensity = inputs.initialDensity;
    const double reynoldsNumber = inputs.reynoldsNumber;
    const double tau = inputs.tau;
    const double maxError = inputs.maxError;

    // Pre-Loop Calculations
    std::cout << "CASE INFO --------\n";
    std::cout << " | Input File: " << inputFileName << "\n";
    std::cout << " | Domain Size: " << numX << " x " << numY << "\n";
    std::cout << " | Maximum Itterations: " << maxItt << "\n";
    std::cout << " | Reynolds Number: " << reynoldsNumber << "\n";
    std::cout << " | Relaxation Constant: " << tau << "\n";
    std::cout << " | Maximum Error: " << maxError << "\n";

    // Solver Variable Definitions
    Matrix2D<double> rho(numY, numX);
    rho.fill(initialDensity); // note to self: make initial fill constructor

    //U

    //f
    //feq
    //fnew

    // Main Loop
    std::cout << "MAIN LOOP --------\n";
    for (int itt = 1; itt <= maxItt; ++itt) {
        if (10 % itt) { std::cout << " | Itt: " << itt << "\n"; }

        // Streaming

        // Collision

        // BGK + Update

    }
    // Output Results


    return 0;
}
