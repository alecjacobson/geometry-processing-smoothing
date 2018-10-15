#include "smooth.h"
#include "cotmatrix.h"
#include "massmatrix.h"
#include <igl/edge_lengths.h>
#include <iostream>

typedef Eigen::Triplet<double> T;
void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
  // Find Laplacian
  Eigen::MatrixXd l;
  igl::edge_lengths(V, F, l);
  Eigen::SparseMatrix<double> L;
  cotmatrix(l, F, L);

  // Find M
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> M;
  massmatrix(l, F, M);
  std::vector<T> triplet_list;
  for(int i = 0; i < M.rows(); i ++) {
    triplet_list.push_back(T(i, i, M.diagonal()(i)));
  }
  Eigen::SparseMatrix<double> sparse_M(M.rows(), M.rows());
  sparse_M.setFromTriplets(triplet_list.begin(), triplet_list.end());

  // Solve system:
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(-lambda * L + sparse_M);
  U = solver.solve(M * G);
}
