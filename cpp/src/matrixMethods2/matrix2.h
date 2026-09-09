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

    int* getDimensions() const;
    double& at(const int& coords);
    void toFile(const std::string& filepath) const;
private:
    int size;
    const int* DIM; // first value is length
    type* data;

    int calcLocation(const int& coords);
    bool equalSizeCheck(const int size1, const int size2);
};

#include "matrix2.tpp"

#endif
