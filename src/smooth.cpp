#include "smooth.h"

#include "massmatrix.h"
#include "cotmatrix.h"
#include "igl/edge_lengths.h"

void smooth(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	const Eigen::MatrixXd & G,
	double lambda,
	Eigen::MatrixXd & U)
{
	Eigen::MatrixXd l;
	igl::edge_lengths(V, F, l);

	Eigen::SparseMatrix<double> L;
	Eigen::DiagonalMatrix<double, Eigen::Dynamic> M;
	cotmatrix(l, F, L);
	massmatrix(l, F, M);
	
	Eigen::SparseMatrix<double> mass;
	{
		std::vector<Eigen::Triplet<double>> entries;
		entries.reserve(M.diagonal().size());
		for (int i = 0; i < M.diagonal().size(); i++) {
			entries.emplace_back(i, i, M.diagonal()(i));
		}
		mass.resize(entries.size(), entries.size());
		mass.setFromTriplets(entries.begin(), entries.end());
	}

	Eigen::SparseMatrix<double> A = (mass - lambda*L);
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> cholesky(A);
	Eigen::MatrixXd b = (M*G).eval();

	U = cholesky.solve(b);
}
