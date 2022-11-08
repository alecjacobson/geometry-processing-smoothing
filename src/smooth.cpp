#include "smooth.h"
#include "cotmatrix.h"
#include "massmatrix.h"
#include "igl/edge_lengths.h"
#include <Eigen/Dense>
#include <iostream>

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
  // Replace with your code
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> M;
  Eigen::MatrixXd l;
  igl::edge_lengths(V, F, l);
  massmatrix(l, F, M);
  Eigen::SparseMatrix<double> L; 
  cotmatrix(l, F, L);
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
  Eigen::MatrixXd Mu = M * G;
  Eigen::SparseMatrix<double> A(M.rows(), M.cols());
  for (int i = 0; i < M.rows();i++) {
  	A.coeffRef(i, i) = M.diagonal()[i];
  }
  A -= lambda * L;
  solver.compute(A);
  U = solver.solve(Mu);
  std::cout << "smooth!" << std::endl;
}
