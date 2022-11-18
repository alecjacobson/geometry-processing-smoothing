#include "massmatrix.h"
#include <igl/doublearea.h>
#include <iostream>

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  // Calculate Area
  Eigen::MatrixXd dblA;
  igl::doublearea(l, dblA);
  dblA = dblA / 2.0;

  // Constuct M:
  int num_V = F.maxCoeff() + 1;
  M.setZero(num_V);
  for(int i = 0; i < F.rows(); i ++){
    double area = dblA(i);
    for(int j = 0; j < 3; j ++) {
      int cur = F(i, j);
      M.diagonal()(cur) += area / 3.0;
    }
  }
}

