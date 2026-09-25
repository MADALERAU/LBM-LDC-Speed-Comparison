#ifndef MATRIX2_H
#define MATRIX2_H

#include <string>

template <class type>
class Matrix {
public:
    Matrix(const int* dim);
    ~Matrix();
    Matrix(const Matrix& source);
    Matrix& operator=(const Matrix& source);
    Matrix(Matrix&& source);
    Matrix& operator=(Matrix&& source);

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

    // The following two methods are vibe coded, will adjust later
    Matrix<type> multiply(const Matrix<type>& rhs) const;
    Matrix<type> multiplyPages(const Matrix<type>& rhs) const;

    // Also vibe coded, will fix later
    Matrix<type> sumOverDim(int dim);

    const int* getDimensions() const;
    type& at(const int location);
    const type& at(const int location) const;
    type& at(const int* coords);
    const type& at(const int* coords) const;

    // Raw contiguous data access. Bypasses calcLocation()'s per-element
    // overhead, needed for tight numerical kernels (e.g. eqmD2Q9) that
    // walk an entire page/field at a time.
    type* getData();
    const type* getData() const;

    void fill(const type value);
    void toFile(const std::string& filepath) const;
    std::string toString() const;
    int calcLocation(const int* coords) const;

    int getSize() const;
private:
    int size;
    int* DIM; // first value is dimension, others are lengths
    type* data;


    void equalSizeCheck(const int size1, const int size2);
};

#include "matrix2.tpp"

#endif
