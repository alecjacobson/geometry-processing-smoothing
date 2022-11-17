#include "massmatrix.h"

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  int num_vertex = F.maxCoeff() + 1;
  Eigen::VectorXd m(num_vertex);
  Eigen::VectorXd A(F.rows());

  for(int i = 0; i < num_vertex; i++){
    double sum = 0;
    for(int t = 0; t < F.rows(); t++){
      if(i == F(t, 0) or i == F(t, 1) or i == F(t, 2)){
        double s = (l(t, 0) + l(t, 1) + l(t, 2))/2;
        sum += sqrt(s * (s - l(t, 0)) * (s - l(t, 1)) * (s - l(t, 2)));
      }
    }
    m(i) = sum;
  }

  M = m.asDiagonal();
}

