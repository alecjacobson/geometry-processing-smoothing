#include "smooth.h"
#include "cotmatrix.h"
#include "massmatrix.h"
#include "igl/edge_lengths.h"
#include "Eigen/SparseCholesky"

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
  // Replace with your code
  
  // Compute edge lengths
  Eigen::MatrixXd l;
  igl::edge_lengths(V, F, l);

  // Construct cot and mass matrices
  Eigen::SparseMatrix<double> L;
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> M;
  cotmatrix(l, F, L_comp); 
  massmatrix(l, F, M);

  // Solve for U as solution of M G = (M - lambda * L) U;
  Eigen::SparseMatrix<double> I(F.maxCoeff() + 1, F.maxCoeff() + 1);
  I.setIdentity();
  Eigen::SparseMatrix<double> L_lambda = I - lambda * M.inverse() * L;
  L_lambda = M * L_lambda;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(L_lambda);
  for (int i = 0; i < G.cols(); i++) {
	U.col(i) = solver.solve(M * G.col(i));
  }
}
