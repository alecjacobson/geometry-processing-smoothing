#include "massmatrix.h"
#include <igl/doublearea.h>
#include <cmath>
#include <vector>
#include <iostream>

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  Eigen::MatrixXd doubleAreas;
  igl::doublearea(l, doubleAreas);
  Eigen::VectorXd coeffs = Eigen::VectorXd::Zero(F.maxCoeff() + 1);
  double thirdOfArea;
  for (int i = 0; i < F.rows(); i++) {
    thirdOfArea = doubleAreas(i) / 6.0;
    for (int j = 0; j < 3; j++) {
      coeffs(F(i, j)) += thirdOfArea;
    }
  }
  M = coeffs.matrix().asDiagonal();
}

