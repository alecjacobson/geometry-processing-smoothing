#include "massmatrix.h"
#include "igl/doublearea.h"

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  // Add your code here

  // Compute area of each face
  Eigen::VectorXd A;
  igl::doublearea(l, A);
  M.resize(F.maxCoeff() + 1);

  // Mass matrix entry (i, i) is 1/3 * (sum of faces containing vertex i)
  for (int i = 0; i <= F.maxCoeff(); i++) {
	double entry = 0;
	for (int j = 0; j < F.rows(); j++) {
		if (F(j, 0) == i || F(j,1) == i || F(j, 2) == i) {
			entry = entry + A(j);
		}
	}
	M.diagonal()[i] = entry / 3;
  }
}

