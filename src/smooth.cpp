#include "smooth.h"
#include "massmatrix.h"
#include "cotmatrix.h"
#include <igl/edge_lengths.h>
void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
	// Replace with your code
	
	//Calculate edge-lengths
	Eigen::MatrixXd l;
	igl::edge_lengths(V, F, l);

	//Get the cotangent and mass information
	Eigen::DiagonalMatrix<double, Eigen::Dynamic> M;
	massmatrix(l, F,  M);
	Eigen::SparseMatrix<double> L;
	cotmatrix(l, F, L);

	//Convert the mass diagonal matrix, into a sparse matrix
	Eigen::SparseMatrix<double> m;	
	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(M.diagonal().size());
	m.resize(V.rows(), V.rows());
	m.setZero();
	for (int i = 0; i < V.rows(); i++) {
		triplets.emplace_back(i, i, M.diagonal()[i]);
	}
	m.setFromTriplets(triplets.begin(), triplets.end());

	//Solve system using cholesky-decomposition
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> cholesky_decomp;
	cholesky_decomp.compute(m - lambda*L);
	U = cholesky_decomp.solve(m * G);
}
