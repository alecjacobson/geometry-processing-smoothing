#include "smooth.h"
#include "cotmatrix.h"
#include "massmatrix.h"
#include <igl/edge_lengths.h>

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
  // Replace with your code
	// compute edge length
	Eigen::MatrixXd l;
	igl::edge_lengths(V, F, l);

	// compute L and M
	Eigen::SparseMatrix<double> L;
	Eigen::DiagonalMatrix<double, Eigen::Dynamic> M;
	cotmatrix(l, F, L);
	massmatrix(l, F, M);

	// solve U
	Eigen::VectorXd diagM = M.diagonal();
	typedef Eigen::Triplet<double> T;
	typedef Eigen::SparseMatrix<double> SpMat;
	std::vector<T> triplets;
	triplets.reserve(diagM.size()*2);
	for (int i = 0; i < diagM.size(); i++) {
		triplets.push_back(T(i, i, diagM(i)));
	}
	SpMat spaM(M.rows(), M.cols());
	spaM.setFromTriplets(triplets.begin(), triplets.end());
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> chol(spaM - lambda * L);
	U = chol.solve(M*G);
}
