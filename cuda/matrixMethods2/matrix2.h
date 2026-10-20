#ifndef MATRIX2_H
#define MATRIX2_H

#include <string>

// Host matrix with an optional device (GPU) mirror.
//
// The host side behaves exactly as before. The device mirror is a plain
// buffer of the same size, created on demand:
//   allocDevice() - cudaMalloc (no copy)
//   toDevice()    - allocate if needed, then host -> device
//   toHost()      - device -> host
// Kernels use getDeviceData(). Copying a Matrix copies host data only; the
// copy has no device buffer until you ask for one.
template <class type>
class Matrix {
public:
    Matrix(const int* dim);
    ~Matrix();
    Matrix(const Matrix& source);
    Matrix& operator=(const Matrix& source);
    Matrix(Matrix&& source);
    Matrix& operator=(Matrix&& source);

    // Host element-wise operators
    Matrix& operator+=(const Matrix& rhs);
    Matrix& operator+=(const type rhs);
    Matrix& operator-=(const Matrix& rhs);
    Matrix& operator-=(const type rhs);
    Matrix& operator*=(const Matrix& rhs);
    Matrix& operator*=(const type rhs);
    Matrix& operator/=(const Matrix& rhs);
    Matrix& operator/=(const type rhs);

    Matrix operator+(const Matrix& rhs) const;
    Matrix operator-(const Matrix& rhs) const;
    Matrix operator*(const Matrix& rhs) const;
    Matrix operator/(const Matrix& rhs) const;

    Matrix<type> generateTranspose() const;

    const int* getDimensions() const;
    type& at(const int location);
    const type& at(const int location) const;
    type& at(const int* coords);
    const type& at(const int* coords) const;

    // Raw contiguous host data access
    type* getData();
    const type* getData() const;

    // ---- Device mirror ----
    void allocDevice();                 // allocate device buffer if missing
    void toDevice();                    // host -> device (allocates if needed)
    void toHost();                      // device -> host
    type* getDeviceData();
    const type* getDeviceData() const;
    // O(1) exchange of device buffers (replaces an f = fNew copy)
    void swapDevice(Matrix& other);

    void fill(const type value);
    void toFile(const std::string& filepath) const;
    std::string toString() const;
    int calcLocation(const int* coords) const;

    int getSize() const;
private:
    int size;
    int* DIM; // first value is dimension, others are lengths
    type* data;
    type* d_data; // device mirror, nullptr until allocDevice()/toDevice()

    void equalSizeCheck(const int size1, const int size2);
};

#include "matrix2.tpp"

#endif
