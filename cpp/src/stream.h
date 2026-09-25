#ifndef STREAM_H
#define STREAM_H

#include "matrixMethods2/matrix2.h"

template <class type>
void stream(
    Matrix<type>& fNew,
    const Matrix<type>& f,
    const Matrix<type>& rho,
    const double uLid,
    const Matrix<int>& zoneId
);

#include "stream.tpp"

#endif
