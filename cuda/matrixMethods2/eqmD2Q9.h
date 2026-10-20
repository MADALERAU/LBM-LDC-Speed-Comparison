#ifndef EQMD2Q9_H
#define EQMD2Q9_H

#include "matrix2.h"

// D2Q9 equilibrium distribution, GPU version of the MATLAB:
//
//   f_eq = rho .* (1 + pagemtimes(ksi,u)/c_s^2 ...
//                    + pagemtimes(ksi,u).^2/(2*c_s^4) ...
//                    - sum(u.^2)/(2*c_s^2)) .* w;
//
// All Matrix arguments must already be on the device (toDevice()).
//   rho : 2D {numY, numX}
//   ksi : 2D {numDir, 2}
//   u   : 3D {numY, numX, 2}   (page 0 = u_x, page 1 = u_y)
//   w   : 2D {numDir, 1}
//   fEq : 3D {numY, numX, numDir}  (output, preallocated on the device)
//
// One thread per grid node loops over the 9 directions. Pages are stored
// contiguously, so thread n touches element n of every page: fully coalesced.
template <class type>
void eqmD2Q9(
    Matrix<type>& fEq,
    const Matrix<type>& rho,
    const Matrix<type>& ksi,
    const Matrix<type>& u,
    const Matrix<type>& w
);

// Fused BGK collision, equivalent to
//   fEq = eqmD2Q9(rho, ksi, u, w);
//   fNew = fNew*(1 - 1/tau) + fEq/tau;
// but the equilibrium is computed in registers and never stored, which
// removes a full write + read of the 9-page fEq array every iteration.
template <class type>
void bgkCollide(
    Matrix<type>& fNew,
    const Matrix<type>& rho,
    const Matrix<type>& ksi,
    const Matrix<type>& u,
    const Matrix<type>& w,
    const type tau
);

#include "eqmD2Q9.tpp"

#endif
