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
    Matrix<type> Matrix<type>::multiply(const Matrix<type>& rhs) const;
    Matrix<type> Matrix<type>::multiplyPages(const Matrix<type>& rhs) const;

    // Also vibe coded, will fix later
    Matrix<type> Matrix<type>::sumOverDim(int dim);

    const int* getDimensions() const;
    type& Matrix<type>::at(const int location);
    type& Matrix<type>::at(const int* coords);
    void fill(const type value);
    void toFile(const std::string& filepath) const;
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
