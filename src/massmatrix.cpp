#include "massmatrix.h"
#include <igl/sort.h>

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  // Add your code here
  M.resize(F.maxCoeff() + 1);
  M.diagonal() = Eigen::VectorXd::Zero(F.maxCoeff() + 1);

  float a, b, c, area = 0;
  Eigen::ArrayXi face(3);
  Eigen::VectorXf unsorted(3);
  Eigen::VectorXf sorted(3);
  Eigen::VectorXi index;

  int v1, v2, v3 = 0;
  int s;

  for (int i = 0; i < F.rows(); i++)
  {
	  v1 = F(i, 0);
	  v2 = F(i, 1);
	  v3 = F(i, 2);

	  unsorted << l(i, 0), l(i, 1), l(i, 2);
	  igl::sort(unsorted, 1, 0, sorted, index);

	  a = sorted[0];
	  b = sorted[1];
	  c = sorted[2];

	  s = (a + b + c) / 2;
	  //area = sqrt(s*(s-a)*(s-b)*(s-c));
	  area = sqrt((a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c))) / 4 ; // stable Heron's formula
	  M.diagonal()(v1) += area;
	  M.diagonal()(v2) += area;
	  M.diagonal()(v3) += area;
  }

  M.diagonal() /= 3;
}

