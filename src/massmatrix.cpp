#include "massmatrix.h"

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  // Add your code here
  int v = F.maxCoeff() + 1;
	M.resize(v);
	M.setZero();
	Eigen::VectorXd temp(v);
	temp.setZero();

	for (int i = 0; i < F.rows(); i++)
	{
		double a = l(i, 0);
		double b = l(i, 1);
		double c = l(i, 2);
		double s = (a + b + c) / 2.0;
		//by helon's formula 
		double area = sqrt(s * (s - a) * (s - b)* (s - c));		
		temp(F(i, 0)) += area / 3.0;
		temp(F(i, 1)) += area / 3.0;
		temp(F(i, 2)) += area / 3.0;
	}

	M = temp.asDiagonal();

}

