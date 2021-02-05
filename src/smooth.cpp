#include "smooth.h"
#include "cotmatrix.h"
#include "massmatrix.h"

#include "igl/edge_lengths.h"

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
  // Replace with your code
  U.resize(G.rows(), G.cols());

  // construct l
  Eigen::MatrixXd l;
  igl::edge_lengths(V, F, l);

  Eigen::SparseMatrix<double> L;
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> M;
  cotmatrix(l, F, L);
  massmatrix(l, F, M);

  // construct A = M - lambda*L
  Eigen::SparseMatrix<double> A = -lambda*L;
  Eigen::VectorXd diag = M.diagonal();
  for (int i = 0; i < diag.rows(); ++i) {
    A.coeffRef(i,i) += diag(i);
  }

  // construct b = M*G
  Eigen::MatrixXd b = M*G;

  // solve
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
  solver.compute(A);
  U = solver.solve(b);
}
