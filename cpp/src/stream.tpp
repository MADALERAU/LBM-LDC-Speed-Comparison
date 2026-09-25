#include <stdexcept>

template <class type>
void stream(Matrix<type>& fNew, const Matrix<type>& f, const Matrix<type>& rho, const type uLid, const Matrix<int>& zoneId) {
    const int* fDim = fNew.getDimensions();
    const int numY = fDim[1];
    const int numX = fDim[2];
    const int size = numY*numX;

    int coords[2];
    int location;

    int c0; int c1; int c2; int c3; int c4; int c5; int c6; int c7; int c8;
    type rhoB;
    type rhoC;
    type addedP;

    type avg;

    for (int j = 0; j < numY; ++j) {
        for (int i = 0; i < numX; ++i) {
            coords[0] = i;
            coords[1] = j;
            location = zoneId.calcLocation(coords);

            c0 = location;
            c1 = location +   size;
            c2 = location + 2*size;
            c3 = location + 3*size;
            c4 = location + 4*size;
            c5 = location + 5*size;
            c6 = location + 6*size;
            c7 = location + 7*size;
            c8 = location + 8*size;

            switch (zoneId.at(location)) {
                case 0: // interior
                    fNew.at(c0) = f.at(c0           ); // 0
                    fNew.at(c1) = f.at(c1        - 1); // 1
                    fNew.at(c2) = f.at(c2 - numX    ); // 2
                    fNew.at(c3) = f.at(c3        + 1); // 3
                    fNew.at(c4) = f.at(c4 + numX    ); // 4
                    fNew.at(c5) = f.at(c5 - numX - 1); // 5
                    fNew.at(c6) = f.at(c6 - numX + 1); // 6
                    fNew.at(c7) = f.at(c7 + numX + 1); // 7
                    fNew.at(c8) = f.at(c8 + numX - 1); // 8
                    break;
                case 1: // right side  , i + 1        out of bounds --> bounceback
                    fNew.at(c0) = f.at(c0           ); // 0
                    fNew.at(c1) = f.at(c1        - 1); // 1
                    fNew.at(c2) = f.at(c2 - numX    ); // 2
                    fNew.at(c4) = f.at(c4 + numX    ); // 4
                    fNew.at(c5) = f.at(c5 - numX - 1); // 5
                    fNew.at(c8) = f.at(c8 + numX - 1); // 8


                    avg = (fNew.at(c4) - fNew.at(c2))/2.0;
                    fNew.at(c3) = fNew.at(c1);       // 3
                    fNew.at(c6) = fNew.at(c8) + avg; // 6
                    fNew.at(c7) = fNew.at(c5) - avg; // 7
                    break;
                case 2: // top side    , j + 1        out of bounds --> velocity
                    fNew.at(c0) = f.at(c0           ); // 0
                    fNew.at(c1) = f.at(c1        - 1); // 1
                    fNew.at(c2) = f.at(c2 - numX    ); // 2
                    fNew.at(c3) = f.at(c3        + 1); // 3
                    fNew.at(c5) = f.at(c5 - numX - 1); // 5
                    fNew.at(c6) = f.at(c6 - numX + 1); // 6

                    avg = (fNew.at(c3) - fNew.at(c1))/2.0;
                    rhoB = fNew.at(c0) + fNew.at(c1) + fNew.at(c3) + 2*(fNew.at(c2) + fNew.at(c5) + fNew.at(c6));
                    addedP = rhoB*uLid/2;
                    fNew.at(c4) = fNew.at(c2);                // 4
                    fNew.at(c7) = fNew.at(c5) - avg - addedP; // 7
                    fNew.at(c8) = fNew.at(c6) + avg + addedP; // 8
                    break;
                case 3: // left side   , i - 1        out of bounds --> bounceback
                    fNew.at(c0) = f.at(c0           ); // 0
                    fNew.at(c2) = f.at(c2 - numX    ); // 2
                    fNew.at(c3) = f.at(c3        + 1); // 3
                    fNew.at(c4) = f.at(c4 + numX    ); // 4
                    fNew.at(c6) = f.at(c6 - numX + 1); // 6
                    fNew.at(c7) = f.at(c7 + numX + 1); // 7

                    avg = (fNew.at(c4) - fNew.at(c2))/2.0;
                    fNew.at(c1) = fNew.at(c3);       // 1
                    fNew.at(c5) = fNew.at(c7) + avg; // 5
                    fNew.at(c8) = fNew.at(c6) - avg; // 8
                    break;
                case 4: // bottom side , j - 1        out of bounds --> bounceback
                    fNew.at(c0) = f.at(c0           ); // 0
                    fNew.at(c1) = f.at(c1        - 1); // 1
                    fNew.at(c3) = f.at(c3        + 1); // 3
                    fNew.at(c4) = f.at(c4 + numX    ); // 4
                    fNew.at(c7) = f.at(c7 + numX + 1); // 7
                    fNew.at(c8) = f.at(c8 + numX - 1); // 8

                    rhoC = (fNew.at(c3) - fNew.at(c1))/2.0;
                    fNew.at(c2) = fNew.at(c4);       // 2
                    fNew.at(c5) = fNew.at(c7) + avg; // 5
                    fNew.at(c6) = fNew.at(c8) - avg; // 6
                    break;
                case 5: // top-right   , i + 1, j + 1 out of bounds
                    fNew.at(c0) = f.at(c0           ); // 0
                    fNew.at(c1) = f.at(c1        - 1); // 1
                    fNew.at(c2) = f.at(c2 - numX    ); // 2
                    fNew.at(c5) = f.at(c5 - numX - 1); // 5

                    rhoC = (rho.at(c0 - 1) + rho.at(c0 - numX))/2.0;
                    fNew.at(c3) = fNew.at(c1); // 3
                    fNew.at(c4) = fNew.at(c2); // 4
                    fNew.at(c7) = fNew.at(c5); // 7

                    fNew.at(c6) = (rhoC - fNew.at(c0) - fNew.at(c1) - fNew.at(c2) - fNew.at(c3) - fNew.at(c4) - fNew.at(c5) - fNew.at(c7))/2.0; // 6
                    fNew.at(c8) = fNew.at(c6); // 8
                    break;
                case 6: // top-left    , i - 1, j + 1 out of bounds
                    fNew.at(c0) = f.at(c0           ); // 0
                    fNew.at(c2) = f.at(c2 - numX    ); // 2
                    fNew.at(c3) = f.at(c3        + 1); // 3
                    fNew.at(c6) = f.at(c6 - numX + 1); // 6

                    rhoC = (rho.at(c0 + 1) + rho.at(c0 - numX))/2.0;
                    fNew.at(c1) = fNew.at(c3); // 1
                    fNew.at(c4) = fNew.at(c2); // 4
                    fNew.at(c8) = fNew.at(c6); // 8

                    fNew.at(c5) = (rhoC - fNew.at(c0) - fNew.at(c1) - fNew.at(c2) - fNew.at(c3) - fNew.at(c4) - fNew.at(c6) - fNew.at(c8))/2.0; // 5
                    fNew.at(c7) = fNew.at(c5); // 7
                    break;
                case 7: // bottom-left , i - 1, j - 1 out of bounds --> bounceback
                    fNew.at(c0) = f.at(c0           ); // 0
                    fNew.at(c3) = f.at(c3        + 1); // 3
                    fNew.at(c4) = f.at(c4 + numX    ); // 4
                    fNew.at(c7) = f.at(c7 + numX + 1); // 7

                    rhoC = (rho.at(c0 + 1) + rho.at(c0 + numX))/2.0;
                    fNew.at(c1) = fNew.at(c3); // 1
                    fNew.at(c2) = fNew.at(c4); // 2
                    fNew.at(c5) = fNew.at(c7); // 5

                    fNew.at(c6) = (rhoC - fNew.at(c0) - fNew.at(c1) - fNew.at(c2) - fNew.at(c3) - fNew.at(c4) - fNew.at(c5) - fNew.at(c7))/2.0; // 6
                    fNew.at(c8) = fNew.at(c6); // 8
                    break;
                case 8: // bottom-right, i + 1, j - 1 out of bounds --> bounceback
                    fNew.at(c0) = f.at(c0           ); // 0
                    fNew.at(c1) = f.at(c1        - 1); // 1
                    fNew.at(c4) = f.at(c4 + numX    ); // 4
                    fNew.at(c8) = f.at(c8 + numX - 1); // 8

                    rhoC = (rho.at(c0 - 1) + rho.at(c0 + numX))/2.0;
                    fNew.at(c2) = fNew.at(c4); // 2
                    fNew.at(c3) = fNew.at(c1); // 3
                    fNew.at(c6) = fNew.at(c8); // 6

                    fNew.at(c5) = (rhoC - fNew.at(c0) - fNew.at(c1) - fNew.at(c2) - fNew.at(c3) - fNew.at(c4) - fNew.at(c6) - fNew.at(c8))/2.0; // 5
                    fNew.at(c7) = fNew.at(c5); // 7
                    break;
                default:
                    break;

            }
        }
    }

    //delete[] fDim; // might cause problems if not copy of pointer
}
