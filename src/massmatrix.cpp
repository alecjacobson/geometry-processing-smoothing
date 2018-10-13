#include "massmatrix.h"

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
    M.resize(F.maxCoeff() + 1);
    M.setZero();

    // Compute area of triangle for i-th face
    auto triangle_area = [&l](int i){
        double e0 = l(i, 0);
        double e1 = l(i, 1);
        double e2 = l(i, 2);
        double s = (e0 + e1 + e2) / 2;
        return sqrt(s*(s-e0)*(s-e1)*(s-e2));
    };

    double A;
    for (int i = 0; i < F.rows(); ++i) {
        A = triangle_area(i);
        for (int j = 0; j < 3; ++j) {
            M.diagonal()[F(i, j)] += A;
        }
    }

    M.diagonal() = M.diagonal() / 3;
}

