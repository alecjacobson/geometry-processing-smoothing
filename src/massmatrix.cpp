#include "massmatrix.h"

using namespace Eigen;

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{	
	int numVert = F.maxCoeff() + 1;
	VectorXd diag(numVert);
	diag.setZero();
	
	auto area = [](auto a, auto b, auto c) 
	{
		double s = 0.5*(a + b + c);
		return sqrt(s*(s - a)*(s - b)*(s - c));
	};

	for (int fIndex = 0; fIndex < F.rows(); ++fIndex) 
	{
		int v0 = F(fIndex, 0), v1 = F(fIndex, 1), v2 = F(fIndex, 2);
		double A = area(l(fIndex, 0), l(fIndex, 1), l(fIndex, 2));
		diag(v0) += A;
		diag(v1) += A;
		diag(v2) += A;
	}

	M = diag.asDiagonal();
}

