#include "massmatrix.h"

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  // Add your code here

  M.resize(F.maxCoeff() + 1);

  float a, b, c, area, area_sum = 0;
  Eigen::ArrayXi face(3);
  Eigen::VectorXf sorted(3);

  for (int i = 0; i < F.rows(); i++)
  {
	  area_sum = 0;
	  for (int j = 0; j < F.rows(); j++)
	  {
		  face = F.row(j);
		  if ((face == j).any())
		  {
			  sorted << l(i, 0), l(i, 1), l(i, 2);
			  std::sort(sorted.data(), sorted.data() + sorted.size());
			  a = sorted[0];
			  b = sorted[1];
			  c = sorted[2];
			  //double s = (a + b + c) / 2; //semiperimeter
			  //area = sqrt(s*(s - a)*(s - b)*(s - c)); //Heron's formula
			  area = 0.25 * sqrt((a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c))); // stable Heron's formula
			  area_sum += area;
		  }
	  }
	  M.diagonal()[i] = (1/3) * area_sum;
  }

}

