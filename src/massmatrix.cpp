#include "massmatrix.h"
#include <ppl.h>

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  const int n = F.maxCoeff();
  const int m = F.rows();
  Concurrency::parallel_for(size_t(0),size_t(n),[m,&l,&F,&M](const int i)
  {
    double sum = 0;
    for (int j = 0; j < m; j++)
      if (i == F(j, 0) || i == F(j, 1) || i == F(j, 2))
      {
        double s = 0.5*(l(j, 0) + l(j, 1) + l(j, 2));
        sum += sqrt(s*(s - l(j, 0))*(s - l(j, 1))*(s - l(j, 2)));
      }
    M.diagonal()(i) = sum / 3;
  }, Concurrency::auto_partitioner());
}