#include "massmatrix.h"
#include <math.h>

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  // Add your code here

  int vertexNum = F.maxCoeff() + 1;
  Eigen::VectorXd areas = Eigen::VectorXd::Zero(vertexNum);
  for (int i = 0;i<F.rows();i++) {
  	double a = l(i, 0);
  	double b = l(i, 1);
  	double c = l(i, 2);
    double s = (a + b + c) / 2;
    double area = sqrt(s * (s - a) * (s - b) * (s - c)) / 3;
    areas(F(i, 0)) += area;
    areas(F(i, 1)) += area;
    areas(F(i, 2)) += area;
  }
  M = areas.asDiagonal();
}

