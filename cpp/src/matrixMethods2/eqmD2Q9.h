#ifndef EQMD2Q9_H
#define EQMD2Q9_H

#include "matrix2.h"

// D2Q9 equilibrium distribution function, equivalent to the MATLAB:
//
//   c_s = 1/sqrt(3);
//   f_eq = rho .* (1 + pagemtimes(ksi,u)/c_s^2 ...
//                    + pagemtimes(ksi,u).^2/(2*c_s^4) ...
//                    - sum(u.^2)/(2*c_s^2)) .* w;
//
// Arguments:
//   rho : 2D matrix, dims {numY, numX}      - density field
//   ksi : 2D matrix, dims {numDir, 2}       - lattice direction vectors
//   u   : 3D matrix, dims {numY, numX, 2}   - velocity field (page 0 = u_x, page 1 = u_y)
//   w   : 2D matrix, dims {numDir, 1}       - lattice weights
//
// Returns:
//   f_eq : 3D matrix, dims {numY, numX, numDir}, page d = the equilibrium
//          distribution for lattice direction d over the whole grid.
//
// Implemented as a single dedicated kernel (not via Matrix's generic
// +/-/* operators) so no intermediate NxN-sized Matrix temporaries are
// allocated: it makes one pass to precompute |u|^2 once, then one
// contiguous, vectorizable pass per lattice direction directly on the
// underlying data buffers - O(numDir * numY * numX) work with O(numY*numX)
// extra memory (the |u|^2 scratch buffer), on top of the output itself.
template <class type>
Matrix<type> eqmD2Q9(
    const Matrix<type>& rho,
    const Matrix<type>& ksi,
    const Matrix<type>& u,
    const Matrix<type>& w
);

#include "eqmD2Q9.tpp"

#endif
