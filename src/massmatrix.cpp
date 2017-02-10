#include "massmatrix.h"

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
	int m = F.rows();
	for (int t = 0; t < m; ++t) {
		Eigen::RowVectorXi triangle = F.row(t);
		int i = triangle(0);
		int j = triangle(1);
		int k = triangle(2);
		double a = l(t, 0);
		double b = l(t, 1);
		double c = l(t, 2);
		double s = 0.5*(a + b + c);
		double area = sqrt(s*(s - a)*(s - b)*(s - c))/3;
		M.diagonal()(i) += area;
		M.diagonal()(j) += area;
		M.diagonal()(k) += area;
	}
	
}

