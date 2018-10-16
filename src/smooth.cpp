#include "smooth.h"
#include "cotmatrix.h"
#include "massmatrix.h"
#include <igl/edge_lengths.h>
#include <iostream>
#include <vector>
#include <Eigen/SparseCholesky>

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
 
  int v = V.rows();
  Eigen::MatrixXd lengths;
  Eigen::SparseMatrix<double> L(v, v), M(v, v);
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> Mdiag(v);
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  
  igl::edge_lengths(V, F, lengths);
  cotmatrix(lengths, F, L);
  massmatrix(lengths, F, Mdiag);

  // need to turn diagonal M into sparse matrix...
  Eigen::VectorXd Mcoeffs = Mdiag.diagonal();
  std::vector<Eigen::Triplet<double>> triplets;
  for (int i = 0; i < v; i++) {
    triplets.push_back(Eigen::Triplet<double>(i, i, Mcoeffs(i)));
  }
  M.setFromTriplets(triplets.begin(), triplets.end());
  
  // solve
  solver.compute(M - lambda * L);
  U = solver.solve(M * G);
}
