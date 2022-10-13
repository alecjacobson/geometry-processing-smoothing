#include "smooth.h"
#include "cotmatrix.h"
#include "massmatrix.h"
#include <igl/edge_lengths.h>
#include <Eigen\src\Cholesky\LLT.h>

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
  // Replace with your code

  Eigen::MatrixXd l;
  igl::edge_lengths(V, F, l);
  
  Eigen::SparseMatrix<double> L;
  cotmatrix(l, F, L);

  Eigen::DiagonalMatrix<double, Eigen::Dynamic> M;
  massmatrix(l, F, M);

  Eigen::SparseMatrix<double> A;
  A = -lambda * L;


  for (int i = 0; i < A.rows(); i++)
  {
	  A.coeffRef(i, i) += M.diagonal()(i);
  }
  //Eigen::MatrixXd A = M - lambda * L;
  Eigen::LLT<Eigen::MatrixXd> llt;
  llt.compute(A);

  Eigen::MatrixXd b = M*G;
  U = llt.solve(b);

   
}
