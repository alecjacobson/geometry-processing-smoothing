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
  // x = U(t+1); b = MU; A = M + lambda*L;
  Eigen::MatrixXd A = lambda * L;
  //Eigen::MatrixXd A = M + (lambda * L);
  Eigen::VectorXd b = M*U;
  U = G;
}
