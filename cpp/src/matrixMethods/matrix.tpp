#include "matrix.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>

// MATRIX
template <class type>
Matrix<type>::Matrix(int size) : SIZE(size), data(new type[size]) {}
// Testing some bullshit
template <class type>
Matrix<type>::Matrix(type* existingData, int size) : data(existingData), SIZE(size) {}
template <class type>
Matrix<type>::~Matrix() { delete[] data; }

template <class type>
int Matrix<type>::getSize() const { return SIZE; }
template <class type>
type Matrix<type>::at(int location) const { return data[location]; }

template <class type>
void Matrix<type>::set(int location, type value) { data[location] = value; }
template <class type>
void Matrix<type>::fill(type value) {
    std::fill(data, data+SIZE, value);
}

// MATRIX2D
template <class type>
Matrix2D<type>::Matrix2D(int numRows, int numCols) :
    NUM_ROWS(numRows),
    NUM_COLS(numCols),
    Matrix<type>(numRows*numCols)
    {}
template <class type>
Matrix2D<type>::Matrix2D(const Matrix2D<type>& source) :
    Matrix<type>(source.getSize()),
    NUM_ROWS(source.getNumRows()),
    NUM_COLS(source.getNumCols()) {
    std::copy(source.data, source.data+source.getSize(), this->data);
}
template <class type>
Matrix2D<type>& Matrix2D<type>::operator=(const Matrix2D<type>& source) {
    if (this != &source) {
        std::copy(source.data, source.data+source.getSize(), data);
    }
    return *this;
}
template <class type>
Matrix2D<type>::Matrix2D(Matrix2D<type>&& source) :
    NUM_ROWS(source.NUM_ROWS),
    NUM_COLS(source.NUM_COLS),
    //data(source.data),
    //SIZE(source.SIZE)
    Matrix<type>(source.data, source.SIZE)
    { source.data = nullptr; }
template <class type>
Matrix2D<type>& Matrix2D<type>::operator=(Matrix2D<type>&& source) {
    delete data;
    data = source.data;
    source.data = nullptr;
    return *this;
}
template <class type>
int Matrix2D<type>::getNumRows() const { return NUM_ROWS; }
template <class type>
int Matrix2D<type>::getNumCols() const { return NUM_COLS; }
template <class type>
type Matrix2D<type>::at(int row, int col) const { return this->data[calcLocation(row, col)]; }
template <class type>
void Matrix2D<type>::set(int row, int col, type value) { this->data[calcLocation(row, col)] = value; }
template <class type>
void Matrix2D<type>::toFile(const std::string& filepath) const {
    std::ofstream out(filepath);
    if (!out) {
        std::cerr << "Failed to open file: " << filepath << "\n";
        return;
    }
    for (int j = 0; j < NUM_ROWS; ++j) {
        for (int i = 0; i < NUM_COLS; ++i) {
            out  << std::setprecision(15) << Matrix::data[calcLocation(j, i)];
            if (i < NUM_COLS - 1)
                out << " ";
        }
        out << "\n";
    }
    out.close();
}
template <class type>
Matrix2D<type> Matrix2D<type>::operator+=(const Matrix2D<type>& rhs) {
    for (int i = 0; i < SIZE; ++i) {
        this->data[i] += rhs.data[i];
    }
    return *this;
}
template <class type>
Matrix2D<type> Matrix2D<type>::operator+=(const type rhs) {
    for (int i = 0; i < SIZE; ++i) {
        this->data[i] += rhs;
    }
    return *this;
}
template <class type>
Matrix2D<type> Matrix2D<type>::operator-=(const Matrix2D<type>& rhs) {
    for (int i = 0; i < SIZE; ++i) {
        this->data[i] -= rhs.data[i];
    }
    return *this;
}
template <class type>
Matrix2D<type> Matrix2D<type>::operator-=(const type rhs) {
    for (int i = 0; i < SIZE; ++i) {
        this->data[i] -= rhs;
    }
    return *this;
}
template <class type>
Matrix2D<type> Matrix2D<type>::operator*=(const Matrix2D<type>& rhs) {
    for (int i = 0; i < SIZE; ++i) {
        this->data[i] *= rhs.data[i];
    }
    return *this;
}
template <class type>
Matrix2D<type> Matrix2D<type>::operator*=(const type rhs) {
    for (int i = 0; i < SIZE; ++i) {
        this->data[i] *= rhs;
    }
    return *this;
}
template <class type>
Matrix2D<type> Matrix2D<type>::operator/=(const Matrix2D<type>& rhs) {
    for (int i = 0; i < SIZE; ++i) {
        this->data[i] /= rhs.data[i];
    }
    return *this;
}
template <class type>
Matrix2D<type> Matrix2D<type>::operator/=(const type rhs) {
    for (int i = 0; i < SIZE; ++i) {
        this->data[i] /= rhs;
    }
    return *this;
}

template <class type>
Matrix2D<type> Matrix2D<type>::operator+(const Matrix2D<type>& rhs) const {
    Matrix2D<type> output(NUM_ROWS, NUM_COLS);
    for (int i = 0; i < SIZE; ++i) {
        output.data[i] = this->data[i] + rhs.data[i];
    }
    return output;
}
template <class type>
Matrix2D<type> Matrix2D<type>::operator-(const Matrix2D<type>& rhs) const {
    Matrix2D<type> output(NUM_ROWS, NUM_COLS);
    for (int i = 0; i < SIZE; ++i) {
        output.data[i] = this->data[i] - rhs.data[i];
    }
    return output;
}
template <class type>
Matrix2D<type> Matrix2D<type>::operator*(const Matrix2D<type>& rhs) const {
    Matrix2D<type> output(NUM_ROWS, NUM_COLS);
    for (int i = 0; i < SIZE; ++i) {
        output.data[i] = this->data[i] * rhs.data[i];
    }
    return output;
}
template <class type>
Matrix2D<type> Matrix2D<type>::operator/(const Matrix2D<type>& rhs) const {
    Matrix2D<type> output(NUM_ROWS, NUM_COLS);
    for (int i = 0; i < SIZE; ++i) {
        output.data[i] = this->data[i] / rhs.data[i];
    }
    return output;
}

template <class type>
Matrix2D<type> Matrix2D<type>::mult(const Matrix2D<type>& rhs) const {
    if (NUM_COLS != rhs.NUM_ROWS) {
        throw std::invalid_argument("Incompatible dimensions for matrix multiplication");
    }

    Matrix2D<type> output(NUM_ROWS, rhs.NUM_COLS);

    for (int j = 0; j < NUM_ROWS; ++j) {
        for (int i = 0; i < rhs.NUM_COLS; ++i) {

            type sum = 0;

            for (int k = 0; k < NUM_COLS; ++k) {
                sum += this->data[calcLocation(j, k)]*rhs.data[rhs.calcLocation(k, i)];
            }

            output.data[output.calcLocation(j, i)] = sum;

        }
    }

    return output;
}
template <class type>
void Matrix2D<type>::mult(const Matrix2D<type>& lhs, const Matrix2D<type>& rhs) {
    if (lhs.NUM_COLS != rhs.NUM_ROWS || lhs.NUM_ROWS != NUM_ROWS || rhs.NUM_COLS != NUM_COLS) {
        throw std::invalid_argument("Incompatible dimensions for matrix multiplication");
    }

    for (int j = 0; j < NUM_ROWS; ++j) {
        for (int i = 0; i < NUM_COLS; ++i) {

            type sum = 0;

            for (int k = 0; k < lhs.NUM_COLS; ++k) {
                sum += lhs.data[lhs.calcLocation(j, k)]*rhs.data[rhs.calcLocation(k, i)];
            }

            this->data[calcLocation(j, i)] = sum;

        }
    }
}

template <class type>
int Matrix2D<type>::calcLocation(int row, int col) const {
    return NUM_COLS*row + col;
    //return NUM_ROWS*col + row;
}

// MATRIX3D
template <class type>
Matrix3D<type>::Matrix3D(int numPages, int numRows, int numCols) : NUM_PAGES(numPages), Matrix2D<type>(numRows, numCols) { SIZE = NUM_PAGES*NUM_ROWS*NUM_COLS; }
template <class type>
int Matrix3D<type>::getNumPages() const { return NUM_PAGES; }
template <class type>
type Matrix3D<type>::at(int page, int row, int col) const {
    return this->data[calcLocation(page, row, col)]
}
template <class type>
void Matrix3D<type>::set(int page, int row, int col, type value) {
    this->data[calcLocation(page, row, col)] = value;
}
template <class type>
void Matrix3D<type>::toFile(const std::string& filepath) const {
    std::ofstream out(filepath);
    if (!out) {
        std::cerr << "Failed to open file: " << filepath << "\n";
        return;
    }

    for (int k = 0; k < NUM_PAGES; ++k) {
        for (int j = 0; j < NUM_ROWS; ++j) {
            for (int i = 0; i < NUM_COLS; ++i) {
                out  << std::setprecision(15) << Matrix::data[calcLocation(k, j, i)];
                if (i < NUM_COLS - 1)
                    out << " ";
            }
            out << "\n";
        }
        out << "\n";
    }

    out.close();
}

template <class type>
int Matrix3D<type>::calcLocation(int page, int row, int col) const {
    return page*NUM_COLS*NUM_ROWS + NUM_COLS*row + col;
}
