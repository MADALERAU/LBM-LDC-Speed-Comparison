#ifndef CUDAUTILS_H
#define CUDAUTILS_H

#include <cuda_runtime.h>

#include <stdexcept>
#include <string>

// Wrap every CUDA runtime call: throws with file/line on failure.
#define CUDA_CHECK(call)                                                      \
    do {                                                                      \
        cudaError_t err_ = (call);                                            \
        if (err_ != cudaSuccess) {                                            \
            throw std::runtime_error(                                         \
                std::string("CUDA error: ") + cudaGetErrorString(err_) +      \
                " at " + __FILE__ + ":" + std::to_string(__LINE__));          \
        }                                                                     \
    } while (0)

// Number of lattice directions (D2Q9). Compile-time so the per-direction
// loops inside the kernels fully unroll.
constexpr int NUM_DIR = 9;

// Threads per block for every 1D "one thread per grid node" kernel.
constexpr int BLOCK_SIZE = 256;

inline int numBlocks(const int n) {
    return (n + BLOCK_SIZE - 1) / BLOCK_SIZE;
}

#endif
