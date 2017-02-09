#include "massmatrix.h"

#include <igl/doublearea.h>
void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
	// Add your code here
	Eigen::VectorXd mass(F.maxCoeff() + 1);
	mass.setZero();
	Eigen::VectorXd area(F.rows());

	igl::doublearea(l, area); // Area of triangle faces in order

	for (int i = 0; i < F.rows(); i++) {
		for (int j = 0; j < 3; j++) {
			mass(F(i, j)) += area(i) /3;
		}
	}

	M = mass.asDiagonal();

}

