#include "massmatrix.h"
#include <igl/doublearea.h>

using namespace Eigen;

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
	int f = F.rows();
	int n = F.maxCoeff() + 1;
	int m = 3 * f / 2;

	VectorXd A(f);
	igl::doublearea(l, A);

	double oneSixth = 1.0 / 6;

	VectorXd mVec(n);
	mVec.setZero();

	for(int i = 0; i < n; ++i)
		for (int j = 0; j < 3; ++j)
			mVec(F(i, j)) += oneSixth*A(i);

	M = mVec.asDiagonal();
}

