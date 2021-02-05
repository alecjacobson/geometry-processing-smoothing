#include "massmatrix.h"

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  // Add your code here

  // number of vertices
  int v_num = F.maxCoeff() + 1;
  
  // set size of L
  Eigen::VectorXd diag(v_num);
  diag.setZero();

  // loop through each face, add area/3 to elements of M
  for (int i = 0; i < F.rows(); ++i) {
    double sum = 0;
    const Eigen::RowVector3i& f = F.row(i);
    const Eigen::RowVector3d& e = l.row(i);

    double area = 0.25*std::sqrt(
      (e(0) + e(1) + e(2))*(-e(0) + e(1) + e(2))*
      (e(0) - e(1) + e(2))*(e(0) + e(1) - e(2))
    );

    diag(f(0)) += area/3.0;
    diag(f(1)) += area/3.0;
    diag(f(2)) += area/3.0;

  } // end loop i

  M = diag.asDiagonal();
}

