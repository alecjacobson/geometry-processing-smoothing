#include "smooth.h"
#include "cotmatrix.h"
#include "massmatrix.h"
#include "igl/edge_lengths.h"

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
  Eigen::MatrixXd l;
  igl::edge_lengths(V, F, l);
  
  Eigen::SparseMatrix<double> L;
  cotmatrix(l, F, L);

  Eigen::DiagonalMatrix<double,Eigen::Dynamic> M;
  massmatrix(l, F, M);

  // M * G = (M - lambda * L) * U 
  Eigen::SparseMatrix<double> A = -lambda * L;
  for(int i = 0; i < M.rows(); i++) {
  	A.coeffRef(i, i) += M.diagonal()[i];
  }

  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver(A);
  U = solver.solve(M * G);
}
