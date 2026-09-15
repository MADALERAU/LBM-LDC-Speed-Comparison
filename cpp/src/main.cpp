#include "inputHandler.h"
#include "matrixMethods2/matrix2.h"

#include <string>
#include <iostream>
#include <cmath>

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
    // general stuff
    double cs = 1.0/sqrt(3);
    double vis = (tau - 0.5)*cs*cs;
    double u_lid = reynoldsNumber*vis/numX;

    // ksi
    const int KSI_DIM[] = {2, 9, 2};
    Matrix<double> ksi(KSI_DIM);

    ksi.at(0) = 0; ksi.at(1) = 0;
    ksi.at(2) = 1; ksi.at(3) = 0;
    ksi.at(4) = 0; ksi.at(5) = 1;
    ksi.at(6) = -1; ksi.at(7) = 0;
    ksi.at(8) = 0; ksi.at(9) = -1;
    ksi.at(10) = 1; ksi.at(11) = 1;
    ksi.at(12) = -1; ksi.at(13) = 1;
    ksi.at(14) = -1; ksi.at(15) = -1;
    ksi.at(16) = 1; ksi.at(17) = -1;

    // w
    const int W_DIM[] = {2, 9, 1};
    Matrix<double> w(W_DIM);

    w.at(0) = 4/9.0;
    w.at(1) = 1/9.0; w.at(2) = 1/9.0; w.at(3) = 1/9.0; w.at(4) = 1/9.0;
    w.at(5) = 1/36.0; w.at(6) = 1/36.0; w.at(7) = 1/36.0; w.at(8) = 1/36.0;

    // rho
    const int RHO_DIM[] = {2, numY, numX};
    Matrix<double> rho(RHO_DIM);
    rho.fill(initialDensity); // note to self: make initial fill constructor

    //U
    const int U_DIM[] = {3, numY, numX, 2};
    Matrix<double> u(U_DIM);
    u.fill(0.0);

    //f

    //feq
    //fnew

    // TESTIN
    ksi.toFile("output/ksi.dat");
    w.toFile("output/w.dat");
    rho.toFile("output/rho.dat");
    u.toFile("output/u.dat");

    return 0;

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

/*
void eqmD2Q9(Matrix<double>& fEq, const Matrix<double>& rho, const Matrix<double>& ksi, const Matrix<double>& u, const Matrix<double>& w, double cs) {
    Matrix<double> middle(fEq.getDimensions());
    Matrix<double> end(fEq.getDimensions());

    fEq = rho;

    middle = ksi.multiplyPages(u);
    middle /= cs*cs;
    middle += 1;
    Matrix<double> middle2 = middle * middle;
    middle2 /= 2*cs*cs*cs*cs;
    middle += middle2;

    Matrix<double> end = u;
    end*= u;
    end.sumOverDim(3);

    fEq *= middle;


    delete[] middle;
    delete[] middle2;
}*/
