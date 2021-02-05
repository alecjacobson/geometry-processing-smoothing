#include "smooth.h"
#include "massmatrix.h"
#include "cotmatrix.h"

#include "igl/edge_lengths.h"
#include <Eigen/IterativeLinearSolvers>

void smooth(const Eigen::MatrixXd & V, const Eigen::MatrixXi & F,
		const Eigen::MatrixXd & G, double lambda, Eigen::MatrixXd & U) {
	// Replace with your code
// first, calculate face edges;
	Eigen::MatrixXd l;
	igl::edge_lengths(V, F, l);

	// second, calculate L
	Eigen::SparseMatrix<double> L;
	cotmatrix(l, F, L);

	// third , M
	Eigen::DiagonalMatrix<double, Eigen::Dynamic> M;
	massmatrix(l, F, M);

	// solve equation
	// next, we combine G with M
	Eigen::SparseMatrix<double> A = -lambda * L;
	for (int i = 0; i < V.rows(); i++) {
		A.coeffRef(i, i) += M.diagonal()(i);
	}
	// last, we solve Ax = B
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
	solver.compute(A);

	Eigen::MatrixXd B = M * G;
	U = solver.solve(B);
}
