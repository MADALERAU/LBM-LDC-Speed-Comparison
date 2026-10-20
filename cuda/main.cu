#include "inputHandler.h"
#include "cudaUtils.h"
#include "matrixMethods2/matrix2.h"
#include "matrixMethods2/eqmD2Q9.h"
#include "macroscopic.h"
#include "stream.h"

#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>

using clock_type = std::chrono::high_resolution_clock;
using duration_type = std::chrono::duration<double, std::milli>;

// GPU launches are asynchronous. With PROFILE on, we synchronize after each
// section so the per-section timers are meaningful (costs a little speed).
// Set to false for production runs.
constexpr bool PROFILE = true;

static void syncIfProfiling() {
    if (PROFILE) { CUDA_CHECK(cudaDeviceSynchronize()); }
}

void populateZoneId(Matrix<int>& zoneId);

int run() {

    std::cout << "MAIN STARTED -----\n";

    std::string inputFileName = "case0";

    InputValues inputs = InputValues(inputFileName);

    // Read Inputs
    const int numX = inputs.numSideNodes;
    const int numY = numX;
    const int maxItt = inputs.maxItt;
    const double initialDensity = inputs.initialDensity;
    const double reynoldsNumber = inputs.reynoldsNumber;
    const double tau = inputs.tau;
    const double maxError = inputs.maxError;

    // Device info
    cudaDeviceProp prop;
    CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));

    // Pre-Loop Calculations

    // Solver Variable Definitions
    // general stuff
    double cs = 1.0/sqrt(3);
    double vis = (tau - 0.5)*cs*cs;
    double uLid = reynoldsNumber*vis/numX;
    double mach = uLid/cs;

    // Every Matrix is set up on the host first, then uploaded with toDevice().

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
    rho.fill(initialDensity);

    // U
    const int U_DIM[] = {3, numY, numX, 2};
    Matrix<double> u(U_DIM);
    u.fill(0.0);

    // zone map
    Matrix<int> zoneId(RHO_DIM);
    populateZoneId(zoneId);

    ksi.toDevice();
    w.toDevice();
    rho.toDevice();
    u.toDevice();
    zoneId.toDevice();

    // f: initial equilibrium, computed on the device
    const int F_DIM[] = {3, numY, numX, 9};
    Matrix<double> f(F_DIM);
    f.allocDevice();
    eqmD2Q9(f, rho, ksi, u, w);
    f.toHost(); // only so the initial state can be written to file below

    // fNew: same size as f, fully overwritten by stream() every iteration
    Matrix<double> fNew = f;
    fNew.toDevice();

    // TESTIN
    ksi.toFile("output/initials/ksi.dat");
    w.toFile("output/initials/w.dat");
    rho.toFile("output/initials/rho.dat");
    u.toFile("output/initials/u.dat");
    f.toFile("output/initials/f.dat");
    fNew.toFile("output/initials/fNew.dat");
    zoneId.toFile("output/initials/zoneId.dat");

    std::cout << "CASE INFO --------\n";
    std::cout << " | Input File: " << inputFileName << "\n";
    std::cout << " | Device: " << prop.name << "\n";
    std::cout << " | Domain Size: " << numX << " x " << numY << "\n";
    std::cout << " | Maximum Itterations: " << maxItt << "\n";
    std::cout << " | Reynolds Number: " << reynoldsNumber << "\n";
    std::cout << " | Relaxation Constant: " << tau << "\n";
    std::cout << " | Maximum Error: " << maxError << "\n";
    std::cout << " | Mach: " << mach << "\n";

    // Main Loop
    std::cout << "MAIN LOOP --------\n";

    double total_loop_time = 0.0;
    double total_stream_time = 0.0;
    double total_collision_time = 0.0;
    double total_fEq_time = 0.0;

    CUDA_CHECK(cudaDeviceSynchronize());
    auto start_time_now = clock_type::now();

    int sampleFrequency = 100;

    for (int itt = 1; itt <= maxItt; ++itt) {
        auto start_time_in_loop = clock_type::now();

        // Streaming
        auto stream_start = clock_type::now();
        stream(fNew, f, rho, uLid, zoneId);
        syncIfProfiling();
        auto stream_end = clock_type::now();

        // Collision (macroscopic moments: rho, u)
        auto collision_start = clock_type::now();
        macroscopic(rho, u, fNew, ksi);
        syncIfProfiling();
        auto collision_end = clock_type::now();

        // fEq + BGK + Update (fused kernel)
        auto fEq_start = clock_type::now();
        bgkCollide(fNew, rho, ksi, u, w, tau);
        syncIfProfiling();
        auto fEq_end = clock_type::now();

        // f = fNew, as an O(1) pointer swap instead of a 9*N copy.
        // (fNew now holds the old f, which stream() overwrites completely.)
        f.swapDevice(fNew);

        auto end_time_in_loop = clock_type::now();

        total_stream_time +=
        duration_type(stream_end - stream_start).count();

        total_collision_time +=
        duration_type(collision_end - collision_start).count();

        total_fEq_time +=
        duration_type(fEq_end - fEq_start).count();

        total_loop_time +=
        duration_type(end_time_in_loop - start_time_in_loop).count();

        if (itt % sampleFrequency == 0) {
            std::cout << " | Itteration: " << itt << "\n";
        }
    }

    CUDA_CHECK(cudaDeviceSynchronize());
    auto end_time_now = clock_type::now();

    double wall_time =
    duration_type(end_time_now - start_time_now).count();

    std::cout << std::fixed << std::setprecision(4);

    std::cout << "\n===== PERFORMANCE RESULTS =====\n";
    std::cout << "Iterations: " << maxItt << "\n\n";

    std::cout << "Total wall time: " << wall_time << " ms\n";
    std::cout << "Average wall time per iteration: "
    << wall_time / maxItt << " ms\n\n";

    if (PROFILE) {
        std::cout << "Section                  Total (ms)   Avg/iteration (ms)\n";
        std::cout << "Streaming:               "
        << total_stream_time << "      "
        << total_stream_time / maxItt << "\n";

        std::cout << "Collision (rho, u):      "
        << total_collision_time << "      "
        << total_collision_time / maxItt << "\n";

        std::cout << "fEq + BGK + Update:      "
        << total_fEq_time << "      "
        << total_fEq_time / maxItt << "\n";

        std::cout << "Sum of measured sections: "
        << total_stream_time + total_collision_time + total_fEq_time
        << " ms\n";
    }

    // Bring results back to the host for output
    rho.toHost();
    u.toHost();
    rho.toFile("output/density/rhoX.dat");
    u.toFile("output/velocity/uX.dat");

    return 0;
}

int main() {
    try {
        return run();
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
}

void populateZoneId(Matrix<int>& zoneId) {
    const int* ZONE_DIM = zoneId.getDimensions();
    int numY = ZONE_DIM[1];
    int numX = ZONE_DIM[2];

    int coords[2];
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
            else if (i == 0)            { value = 3; } // left side
            else if (i == numX - 1)     { value = 1; } // right side

            coords[0] = i;
            coords[1] = j;
            zoneId.at(coords) = value;
        }
    }
}
