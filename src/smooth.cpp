#include "smooth.h"
#include "igl/edge_lengths.h"
#include "cotmatrix.h"
#include "massmatrix.h"

void smooth(
		const Eigen::MatrixXd & V,
		const Eigen::MatrixXi & F,
		const Eigen::MatrixXd & G,
		double lambda,
		Eigen::MatrixXd & U) {

	// Compute matrix of edge lengths.
	Eigen::MatrixXd edgeLengths;
	edgeLengths.resizeLike(F);
	igl::edge_lengths(V, F, edgeLengths);

	// Compute Laplacian and Mass Matrices
	Eigen::SparseMatrix<double> laplacian;
	Eigen::DiagonalMatrix<double, Eigen::Dynamic> mass;
	cotmatrix(edgeLengths, F, laplacian);
	massmatrix(edgeLengths, F, mass);

	// Treat the system given as a linear system we can solve with Ax = b.
	// I construct A weirdly because of the data-types involved.
	Eigen::MatrixXd b = mass * G;
	Eigen::SparseMatrix<double> A = -lambda * laplacian;
	for (int i = 0; i < mass.rows(); i++) {
		A.coeffRef(i, i) += mass.diagonal()[i];
	}

	// Now solve A to obtain G for the next iteration.
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
	solver.compute(A);

	U.resizeLike(G);
	U = solver.solve(b);
	return;
}
