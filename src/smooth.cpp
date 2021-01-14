#include "smooth.h"
#include "cotmatrix.h"
#include "massmatrix.h"
#include <igl/edge_lengths.h>
#include <Eigen/SparseCholesky>
#include <vector>

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
  // Replace with your code
  // U = G;

  Eigen::MatrixXd l;
  igl::edge_lengths(V, F, l);

  Eigen::SparseMatrix<double> L;
  cotmatrix(l, F, L);

  Eigen::DiagonalMatrix<double,Eigen::Dynamic> M;
  massmatrix(l, F, M);

  // create a sparse version of M as DiagonalMatrix is not able to add or type cast to a SparseMatrix
  Eigen::VectorXd M_vec = M.diagonal();
  Eigen::SparseMatrix<double> M_sparse;
  std::vector<Eigen::Triplet<double> > M_sparse_vector;
  for (int k = 0; k<M_vec.size(); k++) {
    M_sparse_vector.push_back(Eigen::Triplet<double>(k, k, M_vec(k)));
  }
  M_sparse.resize(M_vec.size(), M_vec.size());
  M_sparse.setFromTriplets(M_sparse_vector.begin(), M_sparse_vector.end());
  

  U.resize(G.rows(), G.cols());

  // Cholesky decomposition and solving
  // The usual SimplicalLLT is very very slow. Resorting to LDLT version of Cholesky
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;

  solver.compute(M_sparse - lambda * L);
  U = solver.solve(M_sparse * G);

}
