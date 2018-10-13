#include "massmatrix.h"
#include <igl/doublearea.h>

using namespace Eigen;

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  M.resize(F.maxCoeff() + 1);
  MatrixXd doubleArea;
  igl::doublearea(l, doubleArea);
  for (int i = 0; i < F.rows(); i++) {
    for (int j = 0; j < 3; j++) {
      int v = F(i, j);
      double area = 0.5 * doubleArea(i, 0);
      M.diagonal()[v] += area / 3.0;
    }
  }
}
