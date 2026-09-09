#include "inputHandler.h"
#include "matrixMethods2/matrix2.h"

#include <string>
#include <iostream>

int main() {
    // Matrix testing
    // Matrix2D<double> matrix1(10, 15);
    // matrix1.fill(2.1);
    // Matrix2D<double> matrix2(10, 15);
    // matrix2.fill(5.0);
    // matrix2+=matrix1;
    // matrix2+=3.2;
    //
    // Matrix2D<double> matrix3 = matrix1 + matrix2;
    // Matrix2D<double> matrix4(10, 15);
    // matrix4.fill(2.1);
    // matrix4*=matrix3;
    //
    // Matrix2D<double> matrix5 = matrix4;
    // matrix5/=2;
    //
    // Matrix2D<double> matrix6 = matrix5/matrix1;
    //
    // Matrix2D<double> matrix7(3, 2);
    // matrix7.fill(2.1);
    // Matrix2D<double> matrix8(2, 5);
    // matrix8.fill(2.0);
    //
    // Matrix2D<double> matrix9 = matrix7.mult(matrix8);
    // Matrix2D<double> matrix10(matrix7.getNumRows(), matrix8.getNumCols());
    // matrix10.mult(matrix7, matrix8);
    //
    // Matrix3D<double> matrix11(9, 3, 2);
    // matrix11.fill(2.2);
    //
    // matrix1.toFile("output/matrix1.dat");
    // matrix2.toFile("output/matrix2.dat");
    // matrix3.toFile("output/matrix3.dat");
    // matrix4.toFile("output/matrix4.dat");
    // matrix5.toFile("output/matrix5.dat");
    // matrix6.toFile("output/matrix6.dat");
    // matrix7.toFile("output/matrix7.dat");
    // matrix8.toFile("output/matrix8.dat");
    // matrix9.toFile("output/matrix9.dat");
    // matrix10.toFile("output/matrix10.dat");
    // matrix11.toFile("output/matrix11.data");

    // Matrix2 testing
    int* dimensions1 = new int[3];
    dimensions1[0] = 2;
    dimensions1[1] = 20;
    dimensions1[2] = 30;
    Matrix<double> matrix1(dimensions1&);
    Matrix.toFile("/output/matrix1.dat");

}
