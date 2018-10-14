#include "massmatrix.h"
#include <iostream>
#include <igl/doublearea.h>

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  // Add your code here
	int nv = F.maxCoeff() + 1;
	M.resize(F.maxCoeff() + 1);
	M.setZero();
	Eigen::MatrixXd dblA(F.rows(), 1);
	igl::doublearea(l, dblA);

	// Fill in the diagonal entries of M
	for (int f = 0; f < F.rows(); f++) {
		M.diagonal()[F(f, 0)] += (1.0 / 6.0) * dblA(f);
		M.diagonal()[F(f, 1)] += (1.0 / 6.0) * dblA(f);
		M.diagonal()[F(f, 2)] += (1.0 / 6.0) * dblA(f);
	}
}
