#ifndef MATRIX_H
#define MATRIX_H

#include <string>

template <class type>
class Matrix {
public:
    int getSize() const;
    type at(int location) const;

    void set(int location, type value);
    void fill(type value);
protected:
    Matrix(int size);
    Matrix(type* existingData, int size);
    ~Matrix();

    const int SIZE;
    type* data;

    virtual void toFile(const std::string& filepath) const = 0;
};

template <class type>
class Matrix2D : public Matrix<type> {
public:
    Matrix2D(int numRows, int numCols);
    ~Matrix2D() = default;
    Matrix2D(const Matrix2D& source);
    Matrix2D& operator=(const Matrix2D& source);
    Matrix2D(Matrix2D&& source);
    Matrix2D& operator=(Matrix2D&& source);

    int getNumRows() const;
    int getNumCols() const;
    type at(int row, int col) const;

    void set(int row, int col, type value);

    void toFile(const std::string& filepath) const;

    Matrix2D operator+=(const Matrix2D& rhs);
    Matrix2D operator+=(const type rhs);
    Matrix2D operator-=(const Matrix2D& rhs);
    Matrix2D operator-=(const type rhs);
    Matrix2D operator*=(const Matrix2D& rhs);
    Matrix2D operator*=(const type rhs);
    Matrix2D operator/=(const Matrix2D& rhs);
    Matrix2D operator/=(const type rhs);

    Matrix2D operator+(const Matrix2D& rhs) const;
    Matrix2D operator-(const Matrix2D& rhs) const;
    Matrix2D operator*(const Matrix2D& rhs) const;
    Matrix2D operator/(const Matrix2D& rhs) const;

    Matrix2D mult(const Matrix2D& rhs) const;
    void mult(const Matrix2D& lhs, const Matrix2D& rhs);
protected:
    const int NUM_ROWS;
    const int NUM_COLS;
    int calcLocation(int row, int col) const;
};

template <class type>
class Matrix3D : public Matrix2D<type> {
public:
    Matrix3D(int numPages, int numRows, int numCols);

    int getNumPages() const;
    type at(int page, int row, int col) const;

    void set(int page, int row, int col, type value);

    void toFile(const std::string& filepath) const;
protected:
    const int NUM_PAGES;

    int calcLocation(int page, int row, int col) const;
};

#include "matrix.tpp"

#endif
