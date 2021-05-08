#include "smooth.h"
#include "massmatrix.h"
#include "cotmatrix.h"
#include <igl/edge_lengths.h>
#include <Eigen/SparseCholesky>

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
  // Replace with your code
  // U = G;

  // My code
	Eigen::MatrixXd l;
  igl::edge_lengths(V,F,l);

	Eigen::DiagonalMatrix<double,Eigen::Dynamic> M;
	Eigen::SparseMatrix<double> L;
	massmatrix(l,F,M);
	cotmatrix(l,F,L);

  Eigen::SparseMatrix<double> A;
  A = - L*lambda;
  for (int ii = 0; ii < V.rows(); ii++){
  	A.coeffRef(ii,ii) += M.diagonal()[ii];
  }

  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(A);
	Eigen::MatrixXd b;
	b = M * G;
	U = solver.solve(b);
}
