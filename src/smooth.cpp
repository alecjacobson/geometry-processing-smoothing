#include "smooth.h"
#include "igl/edge_lengths.h"
#include "cotmatrix.h"
#include "massmatrix.h"

using namespace Eigen;

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
	// Replace with your code
	//U = G;
	MatrixXd l(F.rows(), 3);
	igl::edge_lengths(V, F, l);
	
	SparseMatrix<double>  L;
	DiagonalMatrix<double, Eigen::Dynamic>  M;
	cotmatrix(l, F, L);
	massmatrix(l, F, M);
	
	L = -lambda * L;
	for (int i = 0; i < M.diagonal().size(); ++i) 
	{
		L.coeffRef(i,i) += M.diagonal()(i);
	}
	SimplicialLDLT<Eigen::SparseMatrix<double>> cholesky(L);
	U = cholesky.solve(M*G);
}
