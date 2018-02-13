#include "massmatrix.h"

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  // Add your code here

  M.resize(F.maxCoeff() + 1);


  for (int i = 0; i < F.rows(); i++)
  {
	  Eigen::ArrayXi face = F.row(i);
	  if ((face == i).any())
	  {
		  float a = l(i, 0);
		  float b = l(i, 1);
		  float c = l(i, 2);
		  double s = (a + b + c) / 2; //semiperimeter
		  double area = sqrt(s*(s - a)*(s - b)*(s - c)); //Heron's formula
	  }

  }

}

