#include "massmatrix.h"
#include <math.h>

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  M.resize(F.maxCoeff() + 1);
  for (int f = 0; f < F.rows(); f++) {
    double s = l.row(f).sum() / 2;
    double f_A = sqrt(s * (s-l(f,0)) * (s-l(f,1)) * (s-l(f,2)));
    for (int v = 0; v < 3; v++) {
      int i = F(f, v);
      M.diagonal()[i] += f_A / 3;
    }
  }
}
