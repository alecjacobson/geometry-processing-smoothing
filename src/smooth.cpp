#include "smooth.h"
#include "massmatrix.h"
#include "cotmatrix.h"

#include <igl/edge_lengths.h>
#include <Eigen/SparseCholesky>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>

using namespace Eigen;

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
	MatrixXd l;
	SparseMatrix<double> L;
	DiagonalMatrix<double, Eigen::Dynamic> M;
	//SparseMatrix<double> M;

	igl::edge_lengths(V, F, l);
	cotmatrix(l, F, L);
	//igl::cotmatrix(V, F, L);
	massmatrix(l, F, M);
	//igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);

	SparseMatrix<double> A = (-lambda*L);
	VectorXd diag = M.diagonal();

	int n = F.maxCoeff() + 1;
	for (int i = 0; i < n; ++i)
		A.coeffRef(i, i) += diag(i);

	SimplicialLDLT<SparseMatrix<double>, Lower> solver(A);
	U = solver.solve(M*G);

}
