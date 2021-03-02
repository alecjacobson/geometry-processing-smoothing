#include "smooth.h"
#include "cotmatrix.h"
#include "massmatrix.h"
#include "igl/edge_lengths.h"
#include <vector>
#include <Eigen/SparseCholesky>

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
  Eigen::MatrixXd l(F.rows(), 3);
  igl::edge_lengths(V, F, l);

  int n = F.maxCoeff() + 1;
  Eigen::SparseMatrix<double> L(n,n);
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> M;
  cotmatrix(l, F, L);
  massmatrix(l, F, M);
 
  // construct the system matrix A = M - lambda * L
  Eigen::SparseMatrix<double> A(n,n);
  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  for (int i = 0; i < M.rows(); i++) {
    tripletList.push_back(T(i, i, M.diagonal()[i]));
  }
  A.setFromTriplets(tripletList.begin(), tripletList.end());
  A -= lambda * L;
  
  Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> chol(A);
  U = chol.solve(M*G);
}
