#include "cudaUtils.h"

#include <stdexcept>

// Same logic as the CPU version, with Matrix::at() replaced by raw device
// pointer indexing. Matrix pages are contiguous (location + k*size), so
// neighbouring threads access neighbouring addresses: coalesced reads/writes.
// Only boundary-zone threads diverge (a thin ring of the domain).
template <class type>
__global__ void streamKernel(
    type* __restrict__ fNew,
    const type* __restrict__ f,
    const type* __restrict__ rho,
    const type uLid,
    const int* __restrict__ zoneId,
    const int numX,
    const int numY
) {
    const int size = numX * numY;
    const int location = blockIdx.x * blockDim.x + threadIdx.x;
    if (location >= size) return;

    const int c0 = location;
    const int c1 = location +   size;
    const int c2 = location + 2*size;
    const int c3 = location + 3*size;
    const int c4 = location + 4*size;
    const int c5 = location + 5*size;
    const int c6 = location + 6*size;
    const int c7 = location + 7*size;
    const int c8 = location + 8*size;

    type rhoB;
    type rhoC;
    type addedP;
    type avg;

    switch (zoneId[location]) {
        case 0: // interior
            fNew[c0] = f[c0           ]; // 0
            fNew[c1] = f[c1        - 1]; // 1
            fNew[c2] = f[c2 - numX    ]; // 2
            fNew[c3] = f[c3        + 1]; // 3
            fNew[c4] = f[c4 + numX    ]; // 4
            fNew[c5] = f[c5 - numX - 1]; // 5
            fNew[c6] = f[c6 - numX + 1]; // 6
            fNew[c7] = f[c7 + numX + 1]; // 7
            fNew[c8] = f[c8 + numX - 1]; // 8
            break;
        case 1: // right side  , i + 1        out of bounds --> bounceback
            fNew[c0] = f[c0           ]; // 0
            fNew[c1] = f[c1        - 1]; // 1
            fNew[c2] = f[c2 - numX    ]; // 2
            fNew[c4] = f[c4 + numX    ]; // 4
            fNew[c5] = f[c5 - numX - 1]; // 5
            fNew[c8] = f[c8 + numX - 1]; // 8

            avg = (fNew[c4] - fNew[c2]) / type(2);
            fNew[c3] = fNew[c1];       // 3
            fNew[c6] = fNew[c8] + avg; // 6
            fNew[c7] = fNew[c5] - avg; // 7
            break;
        case 2: // top side    , j + 1        out of bounds --> velocity
            fNew[c0] = f[c0           ]; // 0
            fNew[c1] = f[c1        - 1]; // 1
            fNew[c2] = f[c2 - numX    ]; // 2
            fNew[c3] = f[c3        + 1]; // 3
            fNew[c5] = f[c5 - numX - 1]; // 5
            fNew[c6] = f[c6 - numX + 1]; // 6

            avg = (fNew[c3] - fNew[c1]) / type(2);
            rhoB = fNew[c0] + fNew[c1] + fNew[c3] + 2*(fNew[c2] + fNew[c5] + fNew[c6]);
            addedP = rhoB*uLid / type(2);
            fNew[c4] = fNew[c2];                // 4
            fNew[c7] = fNew[c5] - avg - addedP; // 7
            fNew[c8] = fNew[c6] + avg + addedP; // 8
            break;
        case 3: // left side   , i - 1        out of bounds --> bounceback
            fNew[c0] = f[c0           ]; // 0
            fNew[c2] = f[c2 - numX    ]; // 2
            fNew[c3] = f[c3        + 1]; // 3
            fNew[c4] = f[c4 + numX    ]; // 4
            fNew[c6] = f[c6 - numX + 1]; // 6
            fNew[c7] = f[c7 + numX + 1]; // 7

            avg = (fNew[c4] - fNew[c2]) / type(2);
            fNew[c1] = fNew[c3];       // 1
            fNew[c5] = fNew[c7] + avg; // 5
            fNew[c8] = fNew[c6] - avg; // 8
            break;
        case 4: // bottom side , j - 1        out of bounds --> bounceback
            fNew[c0] = f[c0           ]; // 0
            fNew[c1] = f[c1        - 1]; // 1
            fNew[c3] = f[c3        + 1]; // 3
            fNew[c4] = f[c4 + numX    ]; // 4
            fNew[c7] = f[c7 + numX + 1]; // 7
            fNew[c8] = f[c8 + numX - 1]; // 8

            avg = (fNew[c3] - fNew[c1]) / type(2);
            fNew[c2] = fNew[c4];       // 2
            fNew[c5] = fNew[c7] + avg; // 5
            fNew[c6] = fNew[c8] - avg; // 6
            break;
        case 5: // top-right   , i + 1, j + 1 out of bounds
            fNew[c0] = f[c0           ]; // 0
            fNew[c1] = f[c1        - 1]; // 1
            fNew[c2] = f[c2 - numX    ]; // 2
            fNew[c5] = f[c5 - numX - 1]; // 5

            rhoC = (rho[c0 - 1] + rho[c0 - numX]) / type(2);
            fNew[c3] = fNew[c1]; // 3
            fNew[c4] = fNew[c2]; // 4
            fNew[c7] = fNew[c5]; // 7

            fNew[c6] = (rhoC - fNew[c0] - fNew[c1] - fNew[c2] - fNew[c3] - fNew[c4] - fNew[c5] - fNew[c7]) / type(2); // 6
            fNew[c8] = fNew[c6]; // 8
            break;
        case 6: // top-left    , i - 1, j + 1 out of bounds
            fNew[c0] = f[c0           ]; // 0
            fNew[c2] = f[c2 - numX    ]; // 2
            fNew[c3] = f[c3        + 1]; // 3
            fNew[c6] = f[c6 - numX + 1]; // 6

            rhoC = (rho[c0 + 1] + rho[c0 - numX]) / type(2);
            fNew[c1] = fNew[c3]; // 1
            fNew[c4] = fNew[c2]; // 4
            fNew[c8] = fNew[c6]; // 8

            fNew[c5] = (rhoC - fNew[c0] - fNew[c1] - fNew[c2] - fNew[c3] - fNew[c4] - fNew[c6] - fNew[c8]) / type(2); // 5
            fNew[c7] = fNew[c5]; // 7
            break;
        case 7: // bottom-left , i - 1, j - 1 out of bounds --> bounceback
            fNew[c0] = f[c0           ]; // 0
            fNew[c3] = f[c3        + 1]; // 3
            fNew[c4] = f[c4 + numX    ]; // 4
            fNew[c7] = f[c7 + numX + 1]; // 7

            rhoC = (rho[c0 + 1] + rho[c0 + numX]) / type(2);
            fNew[c1] = fNew[c3]; // 1
            fNew[c2] = fNew[c4]; // 2
            fNew[c5] = fNew[c7]; // 5

            fNew[c6] = (rhoC - fNew[c0] - fNew[c1] - fNew[c2] - fNew[c3] - fNew[c4] - fNew[c5] - fNew[c7]) / type(2); // 6
            fNew[c8] = fNew[c6]; // 8
            break;
        case 8: // bottom-right, i + 1, j - 1 out of bounds --> bounceback
            fNew[c0] = f[c0           ]; // 0
            fNew[c1] = f[c1        - 1]; // 1
            fNew[c4] = f[c4 + numX    ]; // 4
            fNew[c8] = f[c8 + numX - 1]; // 8

            rhoC = (rho[c0 - 1] + rho[c0 + numX]) / type(2);
            fNew[c2] = fNew[c4]; // 2
            fNew[c3] = fNew[c1]; // 3
            fNew[c6] = fNew[c8]; // 6

            fNew[c5] = (rhoC - fNew[c0] - fNew[c1] - fNew[c2] - fNew[c3] - fNew[c4] - fNew[c6] - fNew[c8]) / type(2); // 5
            fNew[c7] = fNew[c5]; // 7
            break;
        default:
            break;
    }
}

template <class type>
void stream(Matrix<type>& fNew, const Matrix<type>& f, const Matrix<type>& rho, const type uLid, const Matrix<int>& zoneId) {
    const int* fDim = fNew.getDimensions();
    const int numY = fDim[1];
    const int numX = fDim[2];
    const int size = numY * numX;

    if (!fNew.getDeviceData() || !f.getDeviceData() ||
        !rho.getDeviceData()  || !zoneId.getDeviceData()) {
        throw std::runtime_error("stream: all matrices must be on the device");
    }

    streamKernel<type><<<numBlocks(size), BLOCK_SIZE>>>(
        fNew.getDeviceData(), f.getDeviceData(), rho.getDeviceData(),
        uLid, zoneId.getDeviceData(), numX, numY);
    CUDA_CHECK(cudaGetLastError());
}
