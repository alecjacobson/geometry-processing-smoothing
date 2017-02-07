#include "massmatrix.h"

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  // Add your code here

	int num_verts = F.maxCoeff() + 1;
	Eigen::VectorXd masses(num_verts);
	masses.setZero();

	for (int i = 0; i < F.rows(); i++) {
		//Compute the area of this triangle using heron's formula.
		double a = l(i, 0);
		double b = l(i, 1);
		double c = l(i, 2);
		double s = (a + b + c) / 2;
		double A = std::sqrt(s*(s - a)*(s - b)*(s - c))/3.0;

		for (int j = 0; j < 3; j++) {
			masses(F(i, j)) += A;
		}
	}

	

	M = masses.asDiagonal();
	
}

