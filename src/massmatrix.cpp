#include "massmatrix.h"
#include <igl/doublearea.h>

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
    M.resize(F.maxCoeff() + 1);

    Eigen::VectorXd area;
    igl::doublearea(l, area);

    for (int i = 0; i < M.rows(); i++) {
      double sum = 0;

      for (int j = 0; j < F.rows(); j++) {
          if (F(j, 0) == i || F(j,1) == i || F(j, 2) == i) {
              sum += area(j);
          }
      }

      M.diagonal()[i] = sum / 3;
    }
}

