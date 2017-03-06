#include "massmatrix.h"
#include <math.h>

// Return the area of the triangle with sides a, b and c
double get_area(double a, double b, double c) {
  double s = 0.5 * (a + b + c);
  return sqrt(s * (s - a) * (s - b) * (s - c));
}

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  M.setZero();
  int v0, v1, v2;
  double area;
  for (int i = 0; i < F.rows(); i++) {
    v0 = F(i,0);
    v1 = F(i,1);
    v2 = F(i,2);
    area = get_area(l(i,0), l(i,1), l(i,2)) / (double) 3;
    M.diagonal()[v0] += area;
    M.diagonal()[v1] += area;
    M.diagonal()[v2] += area;
  }
}