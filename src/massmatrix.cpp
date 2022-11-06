#include "massmatrix.h"
#include <igl/doublearea.h>

void massmatrix(
	const Eigen::MatrixXd & l,
	const Eigen::MatrixXi & F,
	Eigen::DiagonalMatrix<double, Eigen::Dynamic> & M)
{
	int n = F.maxCoeff()+1;
	M.resize(n);
	Eigen::MatrixXd A;
	igl::doublearea(l, A);
	for (int i = 0; i < F.rows(); i++) {
		for (int j = 0; j < 3; j++) {
			M.diagonal()[F(i, j)] += A(i) / 6.0;
		}
	}
}

