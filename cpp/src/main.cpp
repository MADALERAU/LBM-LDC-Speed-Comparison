#include "inputHandler.h"
#include "matrixMethods2/matrix2.h"
#include "matrixMethods2/eqmD2Q9.tpp"
#include "stream.h"

#include <string>
#include <iostream>
#include <cmath>

void populateZoneId(Matrix<int>& zoneId);

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
    double uLid = reynoldsNumber*vis/numX;

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
    //const int F_DIM[] = {3, numY, numX, 9};
    //Matrix<double> f(F_DIM);
    //f.fill(0.0);
    Matrix<double> f = eqmD2Q9(rho, ksi, u, w);

    //feq
    //Matrix<double> fEq = eqmD2Q9(rho, ksi, u, w);
    Matrix<double> fEq = f;
    //fnew
    // Matrix<double> fNew(F_DIM);
    // fNew.fill(0.0);
    Matrix<double> fNew = f;

    //zone map
    Matrix<int> zoneId(RHO_DIM);
    populateZoneId(zoneId);


    // TESTIN
    ksi.toFile("output/initials/ksi.dat");
    w.toFile("output/initials/w.dat");
    rho.toFile("output/initials/rho.dat");
    u.toFile("output/initials/u.dat");
    f.toFile("output/initials/f.dat");
    fEq.toFile("output/initials/fEq.dat");
    fNew.toFile("output/initials/fNew.dat");
    zoneId.toFile("output/initials/zoneId.dat");

    for( int j = 0; j < numY; ++j) {
        for (int i = 0; i < numX; ++i) {
            int* coords = new int[3];
            coords[0] = i;
            coords[1] = j;

            double value;
            if (j == 0) {
                if (i == 0) {
                    value = 1.0;
                }
                else if (i == numX - 1) {
                    value = 2.0;
                }
                else {
                    value = 3.0;
                }
            }
            else if(j == numY - 1) {
                if (i == 0) {
                    value = 4.0;
                }
                else if (i == numX - 1) {
                    value = 5.0;
                }
                else {
                    value = 6.0;
                }
            }
            else if (i == 0) {
                value = 7.0;
            }
            else if (i == numX - 1) {
                value = 8.0;
            }
            else {
                value = 0.0;
            }

            coords[2] = 0; f.at(coords) = value;
            coords[2] = 1; f.at(coords) = value;
            coords[2] = 2; f.at(coords) = value;
            coords[2] = 3; f.at(coords) = value;
            coords[2] = 4; f.at(coords) = value;
            coords[2] = 5; f.at(coords) = value;
            coords[2] = 6; f.at(coords) = value;
            coords[2] = 7; f.at(coords) = value;
            coords[2] = 8; f.at(coords) = value;

            delete[] coords;
        }
    }
    f.toFile("output/0.dat");

    fNew.fill(-1.0);
    stream(fNew, f, rho, uLid, zoneId);
    std::cout << "this code ran\n";
    fNew.toFile("output/1.dat");

    return 0;

    // Main Loop
    std::cout << "MAIN LOOP --------\n";
    for (int itt = 1; itt <= maxItt; ++itt) {
        if (10 % itt) { std::cout << " | Itt: " << itt << "\n"; }

        // Streaming
        stream(fNew, f, rho, uLid, zoneId);

        // Collision

        // BGK + Update

        fNew.toFile("output/"+std::to_string(itt)+".dat");
        f = fNew;
    }
    // Output Results


    return 0;
}

void populateZoneId(Matrix<int>& zoneId) {
    const int* ZONE_DIM = zoneId.getDimensions();
    int numY = ZONE_DIM[1];
    int numX = ZONE_DIM[2];

    int* coords = new int[2];
    int value;

    for (int j = 0; j < numY; ++j) {
        for (int i = 0; i < numX; ++i) {
            value = 0;

            if (j == 0) { // bottom of domain
                if (i == 0)             { value = 7; } // bottom-left corner
                else if (i == numX - 1) { value = 8; } // bottom-right corner
                else                    { value = 4; } // bottom side
            }
            else if (j == numY - 1) { // top of domain
                if (i == 0)             { value = 6; } // top-left corner
                else if (i == numX - 1) { value = 5; } // top-right corner
                else                    { value = 2; } // top side
            }
            else if (i == 0)            { value = 3; } // right side
            else if (i == numX - 1)     { value = 1; } // left side

            coords[0] = i;
            coords[1] = j;
            zoneId.at(coords) = value;
        }
    }

    delete[] ZONE_DIM;
    delete[] coords;
}
