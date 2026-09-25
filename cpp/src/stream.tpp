#include <stdexcept>

template <class type>
void stream(Matrix<type>& fNew, const Matrix<type>& f, const Matrix<type>& rho, const double uLid, const Matrix<int>& zoneId) {
    const int* fDim = fNew.getDimensions();
    const int numY = fDim[1];
    const int numX = fDim[2];
    const int size = numY*numX;

    int coords[2];
    int location;
    for (int j = 0; j < numY; ++j) {
        for (int i = 0; i < numX; ++i) {
            coords[0] = i;
            coords[1] = j;
            location = zoneId.calcLocation(coords);

            switch (zoneId.at(location)) {
                case 0:
                    fNew.at(location         ) = f.at(location                    ); // 0
                    fNew.at(location +   size) = f.at(location + size   - 1       ); // 1
                    fNew.at(location + 2*size) = f.at(location + 2*size - numX    ); // 2
                    fNew.at(location + 3*size) = f.at(location + 3*size + 1       ); // 3
                    fNew.at(location + 4*size) = f.at(location + 4*size + numX    ); // 4
                    fNew.at(location + 5*size) = f.at(location + 5*size - numX - 1); // 5
                    fNew.at(location + 6*size) = f.at(location + 6*size - numX + 1); // 6
                    fNew.at(location + 7*size) = f.at(location + 7*size + numX + 1); // 7
                    fNew.at(location + 8*size) = f.at(location + 8*size + numX - 1); // 8
                    break;
                case 1:
                    break;
                case 2:
                    break;
                case 3:
                    break;
                case 4:
                    break;
                case 5:
                    break;
                case 6:
                    break;
                case 7:
                    break;
                case 8:
                    break;
                default:
                    break;

            }
        }
    }
}
