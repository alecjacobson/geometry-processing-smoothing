#include "smooth.h"
#include "cotmatrix.h"
#include "massmatrix.h"
#include <igl/edge_lengths.h>
#include <algorithm>

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
  // construct l - the length of the edges
  Eigen::MatrixXd l(F.rows(), 3);
  igl::edge_lengths(V, F, l);
  
  int num_vert = V.rows();
  Eigen::SparseMatrix<double> L(num_vert, num_vert);
  cotmatrix(l, F, L);
  
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> M(num_vert);
  massmatrix(l, F, M);
  
  // Can't perform arithmetic b/w Diag & Sparse, so sparsify M
  Eigen::SparseMatrix<double> sparseM(num_vert, num_vert);
  sparseM.reserve(num_vert);
  sparseM.setZero();
  std::vector<Eigen::Triplet<double> > tripletList;
  for (int i = 0; i < num_vert; i++) {
    tripletList.push_back(Eigen::Triplet<double>(i, i, M.diagonal()[i]));
  }
  sparseM.setFromTriplets(tripletList.begin(), tripletList.end());
  
  L = lambda * L;
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> cg;
  cg.compute(sparseM - L);
  U = cg.solve(sparseM * G);
}
