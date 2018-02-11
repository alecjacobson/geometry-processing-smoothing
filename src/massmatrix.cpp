#include "massmatrix.h"
#include <math.h> 

double faceArea(double a, double b, double c){
	double sinC = sqrt(4*a*a*b*b - pow(a*a+b*b-c*c,2)) / (2*a*b);
	return a*b*sinC / 2;
}

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  // Add your code here
  int numV = F.maxCoeff()+1;
  M.resize(numV);
  M.setZero();

  double A;
  for (int ii = 0; ii < F.rows(); ii++){
  	A = faceArea(l(ii,0), l(ii,1), l(ii,2));
  	M.diagonal()[F(ii,0)] += A / 3;
  	M.diagonal()[F(ii,1)] += A / 3;
  	M.diagonal()[F(ii,2)] += A / 3;
  }
}

