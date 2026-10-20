#include "cudaUtils.h"

#include <stdexcept>

template <class type>
__global__ void macroscopicKernel(
    type* __restrict__ rho,
    type* __restrict__ u,
    const type* __restrict__ f,
    const type* __restrict__ ksi,
    const int gridSize
) {
    const int n = blockIdx.x * blockDim.x + threadIdx.x;
    if (n >= gridSize) return;

    type sumF = 0;
    type momX = 0;
    type momY = 0;

    #pragma unroll
    for (int k = 0; k < NUM_DIR; ++k) {
        const type fk = f[n + k * gridSize];
        sumF += fk;
        momX += ksi[2 * k]     * fk;
        momY += ksi[2 * k + 1] * fk;
    }

    rho[n] = sumF;
    u[n] = momX / sumF;
    u[n + gridSize] = momY / sumF;
}

template <class type>
void macroscopic(
    Matrix<type>& rho,
    Matrix<type>& u,
    const Matrix<type>& f,
    const Matrix<type>& ksi
) {
    const int gridSize = rho.getSize();

    if (f.getSize() != gridSize * NUM_DIR || u.getSize() != gridSize * 2) {
        throw std::invalid_argument("macroscopic: dimension mismatch");
    }
    if (!rho.getDeviceData() || !u.getDeviceData() ||
        !f.getDeviceData()   || !ksi.getDeviceData()) {
        throw std::runtime_error("macroscopic: all matrices must be on the device");
    }

    macroscopicKernel<type><<<numBlocks(gridSize), BLOCK_SIZE>>>(
        rho.getDeviceData(), u.getDeviceData(),
        f.getDeviceData(), ksi.getDeviceData(), gridSize);
    CUDA_CHECK(cudaGetLastError());
}
