#include "cudaUtils.h"

#include <stdexcept>

namespace {

// Equilibrium for one direction at one node.
template <class type>
__device__ __forceinline__ type eqmValue(
    const type rho, const type ux, const type uy, const type uSq,
    const type ksiX, const type ksiY, const type wd
) {
    const type cs2     = type(1) / type(3);                  // c_s^2 = 1/3
    const type invCs2  = type(1) / cs2;
    const type inv2Cs2 = type(1) / (type(2) * cs2);
    const type inv2Cs4 = type(1) / (type(2) * cs2 * cs2);

    const type dot = ksiX * ux + ksiY * uy;
    return rho * wd * (type(1) + dot * invCs2 + dot * dot * inv2Cs4 - uSq * inv2Cs2);
}

template <class type>
__global__ void eqmD2Q9Kernel(
    type* __restrict__ fEq,
    const type* __restrict__ rho,
    const type* __restrict__ u,
    const type* __restrict__ ksi,
    const type* __restrict__ w,
    const int gridSize
) {
    const int n = blockIdx.x * blockDim.x + threadIdx.x;
    if (n >= gridSize) return;

    const type r  = rho[n];
    const type ux = u[n];
    const type uy = u[n + gridSize];
    const type uSq = ux * ux + uy * uy;

    #pragma unroll
    for (int d = 0; d < NUM_DIR; ++d) {
        fEq[n + d * gridSize] =
            eqmValue(r, ux, uy, uSq, ksi[2 * d], ksi[2 * d + 1], w[d]);
    }
}

template <class type>
__global__ void bgkCollideKernel(
    type* __restrict__ fNew,
    const type* __restrict__ rho,
    const type* __restrict__ u,
    const type* __restrict__ ksi,
    const type* __restrict__ w,
    const type tau,
    const int gridSize
) {
    const int n = blockIdx.x * blockDim.x + threadIdx.x;
    if (n >= gridSize) return;

    const type r  = rho[n];
    const type ux = u[n];
    const type uy = u[n + gridSize];
    const type uSq = ux * ux + uy * uy;
    const type omega = type(1) / tau;

    #pragma unroll
    for (int d = 0; d < NUM_DIR; ++d) {
        const type fEq =
            eqmValue(r, ux, uy, uSq, ksi[2 * d], ksi[2 * d + 1], w[d]);
        const int idx = n + d * gridSize;
        fNew[idx] = fNew[idx] * (type(1) - omega) + omega * fEq;
    }
}

// Shared argument validation (same checks as the CPU version)
template <class type>
int checkEqmArgs(
    const Matrix<type>& rho, const Matrix<type>& ksi,
    const Matrix<type>& u, const Matrix<type>& w
) {
    const int* rhoDim = rho.getDimensions();
    const int* ksiDim = ksi.getDimensions();
    const int* uDim   = u.getDimensions();
    const int* wDim   = w.getDimensions();

    if (rhoDim[0] != 2 || ksiDim[0] != 2 || uDim[0] != 3 || wDim[0] != 2) {
        throw std::invalid_argument(
            "eqmD2Q9: expected rho, ksi, w to be 2D and u to be 3D");
    }
    const int numY = rhoDim[1];
    const int numX = rhoDim[2];
    if (uDim[1] != numY || uDim[2] != numX || uDim[3] != 2) {
        throw std::invalid_argument(
            "eqmD2Q9: u must have dims {numY, numX, 2} matching rho");
    }
    if (ksiDim[1] != NUM_DIR || ksiDim[2] != 2) {
        throw std::invalid_argument("eqmD2Q9: ksi must have dims {9, 2}");
    }
    if (wDim[1] != NUM_DIR || wDim[2] != 1) {
        throw std::invalid_argument("eqmD2Q9: w must have dims {9, 1}");
    }
    if (!rho.getDeviceData() || !ksi.getDeviceData() ||
        !u.getDeviceData()   || !w.getDeviceData()) {
        throw std::runtime_error("eqmD2Q9: inputs must be on the device");
    }
    return numY * numX;
}

} // namespace

template <class type>
void eqmD2Q9(
    Matrix<type>& fEq,
    const Matrix<type>& rho,
    const Matrix<type>& ksi,
    const Matrix<type>& u,
    const Matrix<type>& w
) {
    const int gridSize = checkEqmArgs(rho, ksi, u, w);
    if (fEq.getSize() != gridSize * NUM_DIR || !fEq.getDeviceData()) {
        throw std::invalid_argument(
            "eqmD2Q9: fEq must be preallocated on the device, size numY*numX*9");
    }

    eqmD2Q9Kernel<type><<<numBlocks(gridSize), BLOCK_SIZE>>>(
        fEq.getDeviceData(), rho.getDeviceData(), u.getDeviceData(),
        ksi.getDeviceData(), w.getDeviceData(), gridSize);
    CUDA_CHECK(cudaGetLastError());
}

template <class type>
void bgkCollide(
    Matrix<type>& fNew,
    const Matrix<type>& rho,
    const Matrix<type>& ksi,
    const Matrix<type>& u,
    const Matrix<type>& w,
    const type tau
) {
    const int gridSize = checkEqmArgs(rho, ksi, u, w);
    if (fNew.getSize() != gridSize * NUM_DIR || !fNew.getDeviceData()) {
        throw std::invalid_argument(
            "bgkCollide: fNew must be on the device, size numY*numX*9");
    }

    bgkCollideKernel<type><<<numBlocks(gridSize), BLOCK_SIZE>>>(
        fNew.getDeviceData(), rho.getDeviceData(), u.getDeviceData(),
        ksi.getDeviceData(), w.getDeviceData(), tau, gridSize);
    CUDA_CHECK(cudaGetLastError());
}
