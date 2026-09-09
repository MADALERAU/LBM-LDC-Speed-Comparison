#include "matrix2.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>

// Constors, destructors, etc
template <class type>
Matrix<type>::Matrix(const int* dim) : DIM(dim) {
        size = 1;
        for(int i = 1; i <= dim[0]; ++i) { size*=dim[i]; }
        data = new type[size];
};
template <class type>
Matrix<type>::~Matrix() {
    delete[] DIM;
    delete[] data;
};
template <class type>
Matrix<type>::Matrix(const Matrix<type>& source) :
    size(source.size),
    DIM(source.DIM),
    data(new type[size]) {
    std::copy(source.data, source.data+source.size, this->data);
};
template <class type>
Matrix<type>& Matrix<type>::operator=(const Matrix<type>& source) {
    if (this != &source) {
        std::copy(source.data, source.data+source.size, data);
    }
    return *this;
};
template <class type>
Matrix<type>::Matrix(Matrix<type>&& source) :
    size(source.size),
    DIM(source.DIM),
    data(source.data) {
    source.DIM = nullptr;
    source.data = nullptr;
};
template <class type>
Matrix<type>& Matrix<type>::operator=(Matrix<type>&& source) {
    delete data;
    data = source.data;
    source.data = nullptr;
    return *this;
}

// Operators
template <class type>
Matrix<type>& Matrix<type>::operator+=(const Matrix<type>& rhs) {
    equalSizeCheck(size, rhs.size);
    for (int i = 0; i < size; ++i) {
        this->data[i] += rhs.data[i];
    }
    return *this;
}
template <class type>
Matrix<type>& Matrix<type>::operator+=(const type rhs) {
    for (int i = 0; i < size; ++i) {
        this->data[i] += rhs;
    }
    return *this;
}
template <class type>
Matrix<type>& Matrix<type>::operator-=(const Matrix<type>& rhs) {
    equalSizeCheck(size, rhs.size);
    for (int i = 0; i < size; ++i) {
        this->data[i] -= rhs.data[i];
    }
    return *this;
}
template <class type>
Matrix<type>& Matrix<type>::operator-=(const type rhs) {
    for (int i = 0; i < size; ++i) {
        this->data[i] -= rhs;
    }
    return *this;
}
template <class type>
Matrix<type>& Matrix<type>::operator*=(const Matrix<type>& rhs) {
    equalSizeCheck(size, rhs.size);
    for (int i = 0; i < size; ++i) {
        this->data[i] *= rhs.data[i];
    }
    return *this;
}
template <class type>
Matrix<type>& Matrix<type>::operator*=(const type rhs) {
    for (int i = 0; i < size; ++i) {
        this->data[i] *= rhs;
    }
    return *this;
}
template <class type>
Matrix<type>& Matrix<type>::operator/=(const Matrix<type>& rhs) {
    equalSizeCheck(size, rhs.size);
    for (int i = 0; i < size; ++i) {
        this->data[i] /= rhs.data[i];
    }
    return *this;
}
template <class type>
Matrix<type>& Matrix<type>::operator/=(const type rhs) {
    for (int i = 0; i < ; ++i) {
        this->data[i] /= rhs;
    }
    return *this;
}

template <class type>
Matrix<type> Matrix<type>::operator+(const Matrix<type>& rhs) const {
    equalSizeCheck(size, rhs.size);
    Matrix<type> output(DIM);
    for (int i = 0; i < size; ++i) {
        output.data[i] = this->data[i] + rhs.data[i];
    }
    return output;
}
template <class type>
Matrix<type> Matrix<type>::operator-(const Matrix<type>& rhs) const {
    equalSizeCheck(size, rhs.size);
    Matrix<type> output(DIM);
    for (int i = 0; i < size; ++i) {
        output.data[i] = this->data[i] - rhs.data[i];
    }
    return output;
}
template <class type>
Matrix<type> Matrix<type>::operator*(const Matrix<type>& rhs) const {
    equalSizeCheck(size, rhs.size);
    Matrix<type> output(DIM);
    for (int i = 0; i < size; ++i) {
        output.data[i] = this->data[i] * rhs.data[i];
    }
    return output;
}
template <class type>
Matrix<type> Matrix<type>::operator/(const Matrix<type>& rhs) const {
    equalSizeCheck(size, rhs.size);
    Matrix<type> output(DIM);
    for (int i = 0; i < size; ++i) {
        output.data[i] = this->data[i] / rhs.data[i];
    }
    return output;
}


// Other methods
template <class type>
int* Matrix<type>::getDimensions() const {
    return DIM;
}
template <class type>
double& Matrix<type>::at(const int& coords) {
    return;
}
template <class type>
void Matrix<type>::fill(const type value) {
    std::cout << "sizef: " << size << "\n";
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
                if (i < DIM[2] - 1)
                    out << " ";
            }
            out << "\n";
        }

        //out << "file\n";
    }
    else {
        std::cerr << "Invalid dimensions for file output\n";
    }


    delete[] coords;

    out.close();
}

template <class type>
int Matrix<type>::getSize() const {
    return size;
}
// Private methods
template <class type>
int Matrix<type>::calcLocation(const int* coords) const {
    int location;
    if (DIM[0] == 2) {
        location = DIM[1]*coords[1] + coords[0];
    }
    if (DIM[0] == 3) {
        location = DIM[3]*DIM[2]*coords[2] + DIM[2]*coords[1] + coords[0];
    }
    return location;
}
template <class type>
void Matrix<type>::equalSizeCheck(const int size1, const int size2) {
    if (size1 != size2) {
        std::cerr << "Sizes do not match for operation" << "\n";
    }
}
