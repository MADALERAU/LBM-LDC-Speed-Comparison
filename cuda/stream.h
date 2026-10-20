#ifndef STREAM_H
#define STREAM_H

#include "matrixMethods2/matrix2.h"

// Streaming + boundary conditions (GPU).
//
// One thread per grid node. Each thread pulls its 9 populations from the
// neighbouring nodes of f (so every write to fNew is unique: no races) and
// applies the zone-specific boundary closure from zoneId.
//
// fNew, f, rho, zoneId must all be on the device (toDevice()).
template <class type>
void stream(
    Matrix<type>& fNew,
    const Matrix<type>& f,
    const Matrix<type>& rho,
    const type uLid,
    const Matrix<int>& zoneId
);

#include "stream.tpp"

#endif
