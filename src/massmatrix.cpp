#include "massmatrix.h"
#include <igl/doublearea.h>
#include <iostream>

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{

	std::cout << "Constructing mass matrix..." << std::flush;

	M.resize(F.maxCoeff()+1);

  	// To construct this matrix, we can simply loop through all the faces, and for
	// every vertex encountered, we update the corresponding matrix entry.

	Eigen::VectorXd vec (F.maxCoeff()+1);		
	// Initialize to 0.0
	for (int ii = 0; ii < F.maxCoeff()+1; ii++)
		vec(ii) = 0.0;

	Eigen::VectorXd dblA (F.rows());

	// Compute areas*2.0
	igl::doublearea(l, 0.0, dblA);

	for (int ii = 0; ii < F.rows(); ii++)
	{
		vec(F(ii, 0)) += dblA(ii)/2.0/3.0;
		vec(F(ii, 1)) += dblA(ii)/2.0/3.0;
		vec(F(ii, 2)) += dblA(ii)/2.0/3.0;
	}

	// M = vec.asDiagonal();
	M.diagonal() = vec;

	// for (int ii = 0; ii < F.maxCoeff()+1; ii++)
	// 	std::cout << M(ii) << std::endl;

	// for (int ii = 0; ii < F.maxCoeff()+1; ii++)
	// 	M.insert(ii, ii) = vec(ii);

	std::cout << "done." << std::endl;

}

