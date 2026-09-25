#include <stdexcept>
#include <vector>

template <class type>
Matrix<type> eqmD2Q9(
    const Matrix<type>& rho,
    const Matrix<type>& ksi,
    const Matrix<type>& u,
    const Matrix<type>& w
) {
    const int* rhoDim = rho.getDimensions();
    const int* ksiDim = ksi.getDimensions();
    const int* uDim   = u.getDimensions();
    const int* wDim   = w.getDimensions();

    if (rhoDim[0] != 2 || ksiDim[0] != 2 || uDim[0] != 3 || wDim[0] != 2) {
        throw std::invalid_argument(
            "eqmD2Q9: expected rho, ksi, w to be 2D and u to be 3D"
        );
    }

    const int numY   = rhoDim[1];
    const int numX   = rhoDim[2];
    const int numDir = ksiDim[1];

    if (uDim[1] != numY || uDim[2] != numX || uDim[3] != 2) {
        throw std::invalid_argument(
            "eqmD2Q9: u must have dims {numY, numX, 2} matching rho"
        );
    }
    if (ksiDim[2] != 2) {
        throw std::invalid_argument("eqmD2Q9: ksi must have dims {numDir, 2}");
    }
    if (wDim[1] != numDir || wDim[2] != 1) {
        throw std::invalid_argument(
            "eqmD2Q9: w must have dims {numDir, 1} matching ksi"
        );
    }

    const int gridSize = numY * numX;

    const type cs2     = type(1) / type(3);                  // c_s^2 = 1/3
    const type invCs2  = type(1) / cs2;
    const type inv2Cs2 = type(1) / (type(2) * cs2);
    const type inv2Cs4 = type(1) / (type(2) * cs2 * cs2);

    int fEqDim[4] = {3, numY, numX, numDir};
    Matrix<type> fEq(fEqDim);

    const type* rhoData = rho.getData();
    const type* ksiData = ksi.getData();
    const type* uData   = u.getData();
    const type* wData   = w.getData();
    type* fEqData        = fEq.getData();

    const type* ux = uData;              // page 0: u_x over the whole grid
    const type* uy = uData + gridSize;   // page 1: u_y over the whole grid

    // |u|^2 is the same for every direction, so compute it once rather
    // than re-deriving it inside the direction loop below.
    std::vector<type> uSq(gridSize);
    for (int n = 0; n < gridSize; ++n) {
        uSq[n] = ux[n] * ux[n] + uy[n] * uy[n];
    }

    for (int d = 0; d < numDir; ++d) {
        const type ksiX = ksiData[2 * d];
        const type ksiY = ksiData[2 * d + 1];
        const type wd   = wData[d];

        type* fEqPage = fEqData + d * gridSize;

        // Tight, contiguous, branch-free loop - auto-vectorizes well.
        for (int n = 0; n < gridSize; ++n) {
            const type dot = ksiX * ux[n] + ksiY * uy[n];
            fEqPage[n] = rhoData[n] * wd *
                (type(1) + dot * invCs2 + dot * dot * inv2Cs4 - uSq[n] * inv2Cs2);
        }
    }

    return fEq;
}
