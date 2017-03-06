#include "smooth.h"
#include "massmatrix.h"
#include "cotmatrix.h"
#include <igl/edge_lengths.h>

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
  Eigen::MatrixXd l;
  igl::edge_lengths(V,F,l);

  Eigen::SparseMatrix<double> L,A;
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> M;

  cotmatrix(l,F,L);
  massmatrix(l,F,M);


  A = +lambda * L;
  for(int i = 0; i < M.diagonal().size(); ++i) {
      A.coeffRef(i,i) += M.diagonal()(i);
   }

  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> chol(A);
  U = chol.solve(M * G);


}
