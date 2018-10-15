#include "smooth.h"
#include "cotmatrix.h"
#include "massmatrix.h"
#include <igl/edge_lengths.h>
#include <Eigen/SparseCholesky>

void smooth(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	const Eigen::MatrixXd & G,
	double lambda,
	Eigen::MatrixXd & U)
{
	// get L and M
	Eigen::MatrixXd l;
	igl::edge_lengths(V, F, l);
	Eigen::SparseMatrix<double> L;
	cotmatrix(l, F, L);
	Eigen::DiagonalMatrix<double, Eigen::Dynamic> M;
	massmatrix(l, F, M);

	// convert M from DiagonalMatrix to SparseMatrix
		/* Eigen::SparseMatrix<double>(M) not work. */
	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	tripletList.reserve(M.rows());

	for (int i = 0; i < M.rows(); i++) {tripletList.push_back(T(i, i, M.diagonal()[i]));}
	Eigen::SparseMatrix<double> MS;
	MS.resize(M.rows(), M.rows());
	MS.setFromTriplets(tripletList.begin(), tripletList.end());

	// compute U_{t+1}
	U.resizeLike(G);
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> cholesky(MS - lambda * L);
	U = cholesky.solve(M * G);
}
