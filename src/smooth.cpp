#include "smooth.h"
#include "cotmatrix.h"
#include "massmatrix.h"
#include <Eigen/SparseCholesky>
#include <igl/edge_lengths.h>

using namespace Eigen;

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
  // construct cotmatrix and massMatrix
  MatrixXd l;
  igl::edge_lengths(V, F, l);
  SparseMatrix<double> L;
  DiagonalMatrix<double, Dynamic> M;
  cotmatrix(l, F, L);
  massmatrix(l, F, M);

  // solve linear equation: M * G = (M - lam * L) * U
  SparseMatrix<double> right(M.rows(), M.rows());
  for(int i = 0; i < M.rows(); i++) {
  	right.coeffRef(i, i) += M.diagonal()[i];
  }
  right -= lambda * L;
  SimplicialLDLT<SparseMatrix<double>> solver(right);
  MatrixXd left = M * G;
  U = solver.solve(left);
}