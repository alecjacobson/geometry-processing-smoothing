#include "cotmatrix.h"
#include <iostream>
#include <math.h>

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
  // Add your code here
  int vertexNum = F.maxCoeff() + 1;
  for (int i = 0;i<F.rows();i++) {
  	L.coeffRef(F(i, 0), F(i, 1)) += cotangentL(l(i, 0), l(i, 1), l(i, 2));
  	L.coeffRef(F(i, 1), F(i, 2)) += cotangentL(l(i, 1), l(i, 2), l(i, 0));
  	L.coeffRef(F(i, 2), F(i, 0)) += cotangentL(l(i, 2), l(i, 0), l(i, 1));
    L.coeffRef(F(i, 1), F(i, 0)) = L.coeffRef(F(i, 0), F(i, 1));
    L.coeffRef(F(i, 2), F(i, 1)) = L.coeffRef(F(i, 1), F(i, 2));
    L.coeffRef(F(i, 0), F(i, 2)) = L.coeffRef(F(i, 2), F(i, 0));
  }
  L -= L * Eigen::VectorXd::Ones(L.cols()).asDiagonal();
}

double cotangentL(
	double c,
	double a,
	double b) {
	double cos = (a * a + b * b - c * c) / (2 * a * b);
	double sin = sqrt(4 * a * a * b * b - pow(a * a + b * b - c * c, 2)) / (2 * a * b);
	return sin / cos;
}
