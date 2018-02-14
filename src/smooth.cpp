#include "smooth.h"
#include <igl/edge_lengths.h>
#include "cotmatrix.h"
#include "massmatrix.h"
#include <iostream>
#include<Eigen/SparseLU>

#include <igl/cotmatrix.h>

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
  // // Replace with your code
  // U = G;


	// Compute the edge lengths
	Eigen::MatrixXd l;
	igl::edge_lengths(V, F, l);

	// for (int ii = 0; ii < l.rows(); ii++)
	//	std::cout << l(ii,0) << ", " << l(ii,1) << ", " << l(ii,2) << std::endl;

	// Compute the cotangent laplacian matrix
	Eigen::SparseMatrix<double> L;
	cotmatrix(l, F, L);
	//igl::cotmatrix(V, F, L);

	// Compute the mass matrix
	Eigen::DiagonalMatrix<double,Eigen::Dynamic> M;
	massmatrix(l, F, M);

	std::cout << "Smoothing..." << std::flush;

	// Construct the system matrix
	Eigen::SparseMatrix<double> A = -lambda*L;
 	for (int ii = 0; ii < A.rows(); ii++)
 		A.coeffRef(ii, ii) += M.diagonal()(ii);

	// RHS vector
	Eigen::MatrixXd b = M*G;
	//std::cout << L << std::endl;

	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

  	solver.compute(A);
  	U = solver.solve(b);

  	std::cout << "done." << std::endl;

}
