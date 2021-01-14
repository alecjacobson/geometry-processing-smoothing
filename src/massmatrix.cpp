#include "massmatrix.h"

double get_area(double a, double b, double c) {
    double c_half = 1.0 * (a + b + c) / 2;
		double area = sqrt(c_half * (c_half - a) * (c_half - b)* (c_half - c));
	  return area;
}
void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  // Add your code here

   // Add your code here
  int V =  F.maxCoeff()+1;
  M.resize(V);
  Eigen::VectorXd h;
  h.resize(V);

  for (int i=0; i < F.rows(); i++) {
    double l0,l1,l2;
    l2 = l(i,2);
    l0 = l(i,0);
    l1 = l(i,1);
    double area = get_area(l2,l0,l1);
    h(F(i,0)) = 1.0 * area / 3;
    h(F(i,1)) = 1.0 * area / 3;
    h(F(i,2)) = 1.0 * area / 3;
    
  }
  M.diagonal() = h;
}

