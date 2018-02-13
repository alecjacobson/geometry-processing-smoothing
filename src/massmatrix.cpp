#include "massmatrix.h"
#include <iostream>
#include <cmath>

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  int V = F.maxCoeff() + 1;
  int num_faces = F.rows();
  Eigen::VectorXd masses = Eigen::VectorXd::Zero(V);
  for(int f_idx = 0; f_idx < num_faces; f_idx += 1){
    // get vertex indices
    int i = F(f_idx, 0); int j = F(f_idx, 1); int k = F(f_idx, 2);
    // get area of triangle, and add it to masses(index)
    double a = l(f_idx, 0); double b = l(f_idx, 1); double c = l(f_idx, 2);
    double sp = (a + b + c) / 2; // semiperimeter
    double area = sqrt(sp * (sp - a) * (sp - b) * (sp - c));
    masses(i) += area; masses(j) += area; masses(k) += area;
  }
  masses /= 3;
  M = masses.asDiagonal();
}
