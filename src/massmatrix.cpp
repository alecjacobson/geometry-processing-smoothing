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
	
	M.resize(n);

	VectorXd A(f);
	igl::doublearea(l, A);

	double oneSixth = 1.0 / 6;

	VectorXd mVec(n);
	mVec.setZero();

	for (int i = 0; i < f; ++i)
	{
		for (int j = 0; j < 3; ++j)
			mVec(F(i, j)) += oneSixth*std::abs(A(i));
	}

	M = mVec.asDiagonal();
}
