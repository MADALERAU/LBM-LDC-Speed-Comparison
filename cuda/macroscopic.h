#ifndef MACROSCOPIC_H
#define MACROSCOPIC_H

#include "matrixMethods2/matrix2.h"

// Moments of f (GPU), replacing the three CPU steps
//   rho.sumOverDimIntoMat(f, 3);
//   u.multiplyPages3D(ksiTranspose, f);
//   u /= rho;
// in a single pass over f:  rho = sum_k f_k,  u = (sum_k ksi_k f_k) / rho.
//
//   f   : 3D {numY, numX, 9}
//   ksi : 2D {9, 2}
//   rho : 2D {numY, numX}   (output)
//   u   : 3D {numY, numX, 2} (output)
// All must be on the device.
template <class type>
void macroscopic(
    Matrix<type>& rho,
    Matrix<type>& u,
    const Matrix<type>& f,
    const Matrix<type>& ksi
);

#include "macroscopic.tpp"

#endif
