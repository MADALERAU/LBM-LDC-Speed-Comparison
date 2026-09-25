//#include "matrix2.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>

// Constors, destructors, etc
template <class type>
Matrix<type>::Matrix(const int* dim) {

        DIM = new int[dim[0] + 1];
        for(int i = 0; i <= dim[0]; ++i) { DIM[i] = dim[i]; }

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
    DIM(new int[source.DIM[0] + 1]),
    data(new type[size]) {
    std::copy(source.data, source.data+source.size, this->data);
    std::copy(source.DIM, source.DIM+source.DIM[0]+1, this->DIM);
};
template <class type>
Matrix<type>& Matrix<type>::operator=(const Matrix<type>& source) {
    if (this != &source) {

        if (size != source.size) {
            delete[] data;
            data = new type[source.size];
            size = source.size;
        }

        delete[] DIM;
        DIM = new int[source.DIM[0] + 1];

        std::copy(
            source.DIM,
            source.DIM + source.DIM[0] + 1,
            DIM
        );

        std::copy(
            source.data,
            source.data + source.size,
            data
        );
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

    if (this != &source) {
        delete[] DIM;
        delete[] data;

        size = source.size;
        DIM = source.DIM;
        data = source.data;

        source.size = 0;
        source.DIM = nullptr;
        source.data = nullptr;
    }

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
    for (int i = 0; i < size; ++i) {
        this->data[i] /= rhs;
    }
    return *this;
}

template <class type>
Matrix<type> Matrix<type>::operator+(const Matrix<type>& rhs) const {
    // equalSizeCheck(size, rhs.size);
    // Matrix<type> output(DIM);
    // for (int i = 0; i < size; ++i) {
    //     output.data[i] = this->data[i] + rhs.data[i];
    // }
    // return output;
    Matrix<type> output(*this);
    output += rhs;
    return output;
}
template <class type>
Matrix<type> Matrix<type>::operator-(const Matrix<type>& rhs) const {
    Matrix<type> output(*this);
    output += rhs;
    return output;
}
template <class type>
Matrix<type> Matrix<type>::operator*(const Matrix<type>& rhs) const {
    Matrix<type> output(*this);
    output += rhs;
    return output;
}
template <class type>
Matrix<type> Matrix<type>::operator/(const Matrix<type>& rhs) const {
    Matrix<type> output(*this);
    output += rhs;
    return output;
}

// The following two methods are vibe coded matrix multiplication, I will take a closer look at later
template <class type>
Matrix<type> Matrix<type>::multiply(const Matrix<type>& rhs) const
{
    if (DIM[0] != 2 || rhs.DIM[0] != 2) {
        throw std::invalid_argument(
            "Matrix multiplication requires two 2D matrices"
        );
    }

    // A = m x n
    // B = n x p
    // C = m x p

    const int rowsA = DIM[1];
    const int colsA = DIM[2];

    const int rowsB = rhs.DIM[1];
    const int colsB = rhs.DIM[2];

    if (colsA != rowsB) {
        throw std::invalid_argument(
            "Incompatible matrix dimensions for multiplication"
        );
    }

    int dimensions[3] = {
        2,
        rowsA,
        colsB
    };

    Matrix<type> output(dimensions);

    for (int i = 0; i < rowsA; ++i) {

        for (int j = 0; j < colsB; ++j) {

            type sum = 0;

            for (int k = 0; k < colsA; ++k) {
                sum +=
                this->data[i * colsA + k] *
                rhs.data[k * colsB + j];
            }

            output.data[i * colsB + j] = sum;
        }
    }

    return output;
}
template <class type>
Matrix<type> Matrix<type>::multiplyPages(const Matrix<type>& rhs) const
{
    // A must be 3D
    if (DIM[0] != 3) {
        throw std::invalid_argument(
            "Left-hand matrix must be 3D for multiplyPages()"
        );
    }

    // B must be 2D
    if (rhs.DIM[0] != 2) {
        throw std::invalid_argument(
            "Right-hand matrix must be 2D for multiplyPages()"
        );
    }

    // A = rowsA × inner × pages
    // B = inner × colsB
    //
    // C = rowsA × colsB × pages

    const int rowsA  = DIM[1];
    const int inner  = DIM[2];
    const int pages  = DIM[3];

    const int rowsB  = rhs.DIM[1];
    const int colsB  = rhs.DIM[2];

    // Matrix multiplication requires:
    //
    // A columns == B rows
    //
    if (inner != rowsB) {
        throw std::invalid_argument(
            "Incompatible matrix dimensions for page multiplication"
        );
    }

    int dimensions[4] = {
        3,
        rowsA,
        colsB,
        pages
    };

    Matrix<type> output(dimensions);

    for (int page = 0; page < pages; ++page) {

        for (int i = 0; i < rowsA; ++i) {

            for (int j = 0; j < colsB; ++j) {

                type sum = 0;

                for (int k = 0; k < inner; ++k) {

                    // A[i][k][page]
                    const int indexA =
                    page * rowsA * inner
                    + i * inner
                    + k;

                    // B[k][j]
                    const int indexB =
                    k * colsB
                    + j;

                    sum +=
                    this->data[indexA] *
                    rhs.data[indexB];
                }

                // C[i][j][page]
                const int indexOutput =
                page * rowsA * colsB
                + i * colsB
                + j;

                output.data[indexOutput] = sum;
            }
        }
    }

    return output;
}

// also vibe coded
template <class type>
Matrix<type> Matrix<type>::sumOverDim(int dim) {
    if (dim < 0 || dim >= DIM[0]) {
        throw std::out_of_range("Invalid dimension");
    }

    const int numDims = DIM[0];

    // Cannot squeeze a 1D matrix down to 0 dimensions
    if (numDims == 1) {
        throw std::invalid_argument(
            "Cannot sum over the only dimension of a 1D matrix"
        );
    }

    // Create dimensions for output.
    // The summed dimension is removed.
    int* outputDims = new int[numDims];

    outputDims[0] = numDims - 1;

    int outputDim = 1;

    for (int i = 0; i < numDims; ++i) {

        if (i != dim) {
            outputDims[outputDim] = DIM[i + 1];
            ++outputDim;
        }
    }

    Matrix<type> output(outputDims);

    delete[] outputDims;

    // Number of elements before the summed dimension
    int before = 1;

    for (int i = 0; i < dim; ++i) {
        before *= DIM[i + 1];
    }

    // Size of dimension being summed
    const int sumLength = DIM[dim + 1];

    // Number of elements after the summed dimension
    int after = 1;

    for (int i = dim + 1; i < numDims; ++i) {
        after *= DIM[i + 1];
    }

    /*
     D ata layou*t:

     [before][sum][after]

     For every combination of "before" and "after",
     sum across the "sum" dimension.
     */

    for (int b = 0; b < before; ++b) {

        for (int a = 0; a < after; ++a) {

            type sum = 0;

            for (int s = 0; s < sumLength; ++s) {

                const int index =
                b * sumLength * after
                + s * after
                + a;

                sum += data[index];
            }

            const int outputIndex =
            b * after + a;

            output.data[outputIndex] = sum;
        }
    }

    return output;
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
type& Matrix<type>::at(const int* coords) {
    int location = calcLocation(coords);
    return this->data[location];
}
template <class type>
void Matrix<type>::fill(const type value) {
    //std::cout << "sizef: " << size << "\n";
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
    else if (DIM[0] == 3) {
        for (int k = 0; k < DIM[3]; ++k) {
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
int Matrix<type>::getSize() const {
    return size;
}
// Private methods
template <class type>
int Matrix<type>::calcLocation(const int* coords) const {
    int location;
    // if (DIM[0] == 2) {
    //     location = DIM[1]*coords[1] + coords[0];
    // }
    // if (DIM[0] == 3) {
    //     location = DIM[3]*DIM[2]*coords[2] + DIM[2]*coords[1] + coords[0];
    // }
    // return location;

    if (DIM[0] == 2) {
        return DIM[2] * coords[1] + coords[0];
    }

    if (DIM[0] == 3) {
        return DIM[2] * DIM[3] * coords[2]
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
