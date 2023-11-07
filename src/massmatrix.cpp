#include "massmatrix.h"
#include <igl/doublearea.h>

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  // Add your code here
  int V_num = F.maxCoeff() + 1;
  M.resize(V_num);

  Eigen::VectorXd diag(V_num);

  for (int i = 0; i < F.rows(); i++) {
    double edge0 = l(i, 0);
    double edge1 = l(i, 1);
    double edge2 = l(i, 2);
    double area = 0.25 * sqrt((edge0 + edge1 + edge2) * (-edge0 + edge1 + edge2) * (edge0 - edge1 + edge2) * (edge0 + edge1 - edge2));

    int v0 = F(i, 0), v1 = F(i, 1), v2 = F(i, 2);

    diag(v0) += area / 3;
    diag(v1) += area / 3;
    diag(v2) += area / 3;
  }

  M.diagonal() = diag;

}

