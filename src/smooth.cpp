#include "smooth.h"
#include "cotmatrix.h"
#include "massmatrix.h"
#include <igl/edge_lengths.h>
#include <Eigen/Cholesky>

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
  // Get l - edge lengths
  Eigen::MatrixXd l;
  igl::edge_lengths(V, F, l);
  // Get L - cotangent laplacian
  Eigen::SparseMatrix<double> L;
  cotmatrix(l, F, L);
  // Get M - mass Matrix
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> M;
  massmatrix(l, F, M);
  // Update by solving linear system Ax = b
  // x = U(t+1); b = MU(t); A = M - lambda*L;
  Eigen::SparseMatrix<double> A = -(lambda * L);
  Eigen::VectorXd md = M.diagonal();
  for(int i=0; i < M.rows(); i++){
    A.coeffRef(i,i) += md(i);
  }
  Eigen::VectorXd b = M*G;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(A);
  U = solver.solve(b);
}
