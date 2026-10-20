#include "cudaUtils.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>

// Constors, destructors, etc
template <class type>
Matrix<type>::Matrix(const int* dim) : d_data(nullptr) {

    DIM = new int[dim[0] + 1];
    for(int i = 0; i <= dim[0]; ++i) { DIM[i] = dim[i]; }

    size = 1;
    for(int i = 1; i <= dim[0]; ++i) { size*=dim[i]; }
    data = new type[size];
}
template <class type>
Matrix<type>::~Matrix() {
    delete[] DIM;
    delete[] data;
    cudaFree(d_data); // no-op for nullptr; not checked, destructors must not throw
}
template <class type>
Matrix<type>::Matrix(const Matrix<type>& source) :
    size(source.size),
    DIM(new int[source.DIM[0] + 1]),
    data(new type[size]),
    d_data(nullptr) {
    std::copy(source.data, source.data+source.size, this->data);
    std::copy(source.DIM, source.DIM+source.DIM[0]+1, this->DIM);
}
template <class type>
Matrix<type>& Matrix<type>::operator=(const Matrix<type>& source) {
    if (this != &source) {

        if (size != source.size) {
            delete[] data;
            data = new type[source.size];
            size = source.size;

            // Old device buffer has the wrong size now
            cudaFree(d_data);
            d_data = nullptr;
        }

        delete[] DIM;
        DIM = new int[source.DIM[0] + 1];

        std::copy(source.DIM, source.DIM + source.DIM[0] + 1, DIM);
        std::copy(source.data, source.data + source.size, data);
    }
    return *this;
}
template <class type>
Matrix<type>::Matrix(Matrix<type>&& source) :
    size(source.size),
    DIM(source.DIM),
    data(source.data),
    d_data(source.d_data) {
    source.DIM = nullptr;
    source.data = nullptr;
    source.d_data = nullptr;
}
template <class type>
Matrix<type>& Matrix<type>::operator=(Matrix<type>&& source) {

    if (this != &source) {
        delete[] DIM;
        delete[] data;
        cudaFree(d_data);

        size = source.size;
        DIM = source.DIM;
        data = source.data;
        d_data = source.d_data;

        source.size = 0;
        source.DIM = nullptr;
        source.data = nullptr;
        source.d_data = nullptr;
    }

    return *this;
}

// Operators (host)
template <class type>
Matrix<type>& Matrix<type>::operator+=(const Matrix<type>& rhs) {
    equalSizeCheck(size, rhs.size);
    for (int i = 0; i < size; ++i) { this->data[i] += rhs.data[i]; }
    return *this;
}
template <class type>
Matrix<type>& Matrix<type>::operator+=(const type rhs) {
    for (int i = 0; i < size; ++i) { this->data[i] += rhs; }
    return *this;
}
template <class type>
Matrix<type>& Matrix<type>::operator-=(const Matrix<type>& rhs) {
    equalSizeCheck(size, rhs.size);
    for (int i = 0; i < size; ++i) { this->data[i] -= rhs.data[i]; }
    return *this;
}
template <class type>
Matrix<type>& Matrix<type>::operator-=(const type rhs) {
    for (int i = 0; i < size; ++i) { this->data[i] -= rhs; }
    return *this;
}
template <class type>
Matrix<type>& Matrix<type>::operator*=(const Matrix<type>& rhs) {
    equalSizeCheck(size, rhs.size);
    for (int i = 0; i < size; ++i) { this->data[i] *= rhs.data[i]; }
    return *this;
}
template <class type>
Matrix<type>& Matrix<type>::operator*=(const type rhs) {
    for (int i = 0; i < size; ++i) { this->data[i] *= rhs; }
    return *this;
}
template <class type>
Matrix<type>& Matrix<type>::operator/=(const Matrix<type>& rhs) {
    equalSizeCheck(size, rhs.size);
    for (int i = 0; i < size; ++i) { this->data[i] /= rhs.data[i]; }
    return *this;
}
template <class type>
Matrix<type>& Matrix<type>::operator/=(const type rhs) {
    for (int i = 0; i < size; ++i) { this->data[i] /= rhs; }
    return *this;
}

template <class type>
Matrix<type> Matrix<type>::operator+(const Matrix<type>& rhs) const {
    Matrix<type> output(*this);
    output += rhs;
    return output;
}
template <class type>
Matrix<type> Matrix<type>::operator-(const Matrix<type>& rhs) const {
    Matrix<type> output(*this);
    output -= rhs;
    return output;
}
template <class type>
Matrix<type> Matrix<type>::operator*(const Matrix<type>& rhs) const {
    Matrix<type> output(*this);
    output *= rhs;
    return output;
}
template <class type>
Matrix<type> Matrix<type>::operator/(const Matrix<type>& rhs) const {
    Matrix<type> output(*this);
    output /= rhs;
    return output;
}

template <class type>
Matrix<type> Matrix<type>::generateTranspose() const {

    if (DIM[0] != 2) {
        throw std::invalid_argument("Transpose only accepts 2D matricies");
    }

    const int DIM_T[] = {2, DIM[2], DIM[1]};
    Matrix<type> transpose(DIM_T);

    for (int j = 0; j < DIM[1]; ++j) {
        for (int i = 0; i < DIM[2]; ++i) {
            transpose.data[i * DIM[1] + j] =
            this->data[j * DIM[2] + i];
        }
    }

    return transpose;
}

// Other methods
template <class type>
const int* Matrix<type>::getDimensions() const {
    return DIM;
}
template <class type>
type& Matrix<type>::at(const int location) {
    return this->data[location];
}
template <class type>
const type& Matrix<type>::at(const int location) const {
    return this->data[location];
}
template <class type>
type& Matrix<type>::at(const int* coords) {
    int location = calcLocation(coords);
    return this->data[location];
}
template <class type>
const type& Matrix<type>::at(const int* coords) const {
    int location = calcLocation(coords);
    return this->data[location];
}
template <class type>
type* Matrix<type>::getData() {
    return data;
}
template <class type>
const type* Matrix<type>::getData() const {
    return data;
}

// Device mirror
template <class type>
void Matrix<type>::allocDevice() {
    if (d_data == nullptr) {
        CUDA_CHECK(cudaMalloc(&d_data, size * sizeof(type)));
    }
}
template <class type>
void Matrix<type>::toDevice() {
    allocDevice();
    CUDA_CHECK(cudaMemcpy(d_data, data, size * sizeof(type), cudaMemcpyHostToDevice));
}
template <class type>
void Matrix<type>::toHost() {
    if (d_data == nullptr) {
        throw std::runtime_error("Matrix::toHost called with no device buffer");
    }
    CUDA_CHECK(cudaMemcpy(data, d_data, size * sizeof(type), cudaMemcpyDeviceToHost));
}
template <class type>
type* Matrix<type>::getDeviceData() {
    return d_data;
}
template <class type>
const type* Matrix<type>::getDeviceData() const {
    return d_data;
}
template <class type>
void Matrix<type>::swapDevice(Matrix<type>& other) {
    if (size != other.size) {
        throw std::invalid_argument("swapDevice requires equal sizes");
    }
    std::swap(d_data, other.d_data);
}

template <class type>
void Matrix<type>::fill(const type value) {
    for (int i = 0; i < size; ++i) {
        this->data[i] = value;
    }
}
template <class type>
void Matrix<type>::toFile(const std::string& filepath) const {
    std::ofstream out(filepath);
    if (!out) {
        std::cerr << "Failed to open file: " << filepath << "\n";
        return;
    }

    int* coords = new int[DIM[0]];

    if (DIM[0] == 2) {
        for (int j = 0; j < DIM[1]; ++j) {
            for (int i = 0; i < DIM[2]; ++i) {
                coords[1] = j;
                coords[0] = i;
                out  << std::setprecision(15) << this->data[calcLocation(coords)];
                if (i < DIM[2] - 1) {
                    out << " ";
                }
            }
            out << "\n";
        }
    }
    else if (DIM[0] == 3) {
        for (int k = 0; k < DIM[3]; ++k) {
            for (int j = 0; j < DIM[1]; ++j) {
                for (int i = 0; i < DIM[2]; ++i) {
                    coords[2] = k;
                    coords[1] = j;
                    coords[0] = i;
                    out  << std::setprecision(15) << this->data[calcLocation(coords)];
                    if (i < DIM[2] - 1) {
                        out << " ";
                    }
                }
                out << "\n";
            }
            out << "\n";
        }
    }
    else {
        std::cerr << "Invalid dimensions for file output\n";
    }

    delete[] coords;

    out.close();
}

template <class type>
std::string Matrix<type>::toString() const {
    std::string out = "MATRIX\n";
    out += "   | Size: " + std::to_string(size) + "\n";
    out += "   | Dim: ";
    for (int i = 1; i <= DIM[0]; ++i) {
        out += std::to_string(DIM[i]) + " ";
    }
    out += "\n";

    return out;
}

template <class type>
int Matrix<type>::getSize() const {
    return size;
}
// Private methods
template <class type>
int Matrix<type>::calcLocation(const int* coords) const {
    if (DIM[0] == 2) {
        return DIM[2] * coords[1] + coords[0];
    }

    if (DIM[0] == 3) {
        // Page-major (SoA) layout: page coords[2] is the slowest-varying
        // index, so each page is a contiguous DIM[1]*DIM[2] block. This is
        // also the layout the GPU kernels rely on for coalesced access.
        return DIM[1] * DIM[2] * coords[2]
        + DIM[2] * coords[1]
        + coords[0];
    }

    throw std::invalid_argument("Unsupported number of dimensions");
}
template <class type>
void Matrix<type>::equalSizeCheck(const int size1, const int size2) {
    if (size1 != size2) {
        throw std::invalid_argument(
            "Matrix sizes do not match"
        );
    }
}
