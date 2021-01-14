#include "smooth.h"
#include "cotmatrix.h"
#include "massmatrix.h"
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <iostream>
#include <igl/edge_lengths.h>

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
  // Replace with your code
  // U = G;
  Eigen::SparseMatrix<double> A;
  Eigen::SparseMatrix<double> L;
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> M;
  Eigen::MatrixXd l,B;
  igl::edge_lengths(V,F,l);
  massmatrix(l,F,M);
  cotmatrix(l,F,L);
  A = -lambda * L;
  Eigen::VectorXd h;
  h = M.diagonal();
  for (int i = 0; i < h.rows(); i++) {
    A.coeffRef(i,i) += h(i);
  }
  // B = Mu
  B = M * G;
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::UpLoType::Lower | Eigen::UpLoType::Upper> solver;
  U = solver.compute(A).solve(B);
  }
