#include "smooth.h"
#include "cotmatrix.h"
#include "massmatrix.h"
#include <igl/edge_lengths.h>
#include <Eigen/IterativeLinearSolvers>
#include <iostream> 

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
  // Replace with your code
  Eigen::MatrixXd l;
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> M;
  Eigen::SparseMatrix<double> L;
  igl::edge_lengths(V, F, l);
  cotmatrix(l, F, L);
  massmatrix(l, F, M);
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
  L *= lambda;
  Eigen::SparseMatrix<double> M_1 = -L;

  for (int i = 0; i< M.rows(); i++){
    M_1.coeffRef(i,i) += M.diagonal()(i);
  }

  Eigen::MatrixXd M_2 = M*G;
  solver.compute(M_1);
  U = solver.solve(M_2);

}