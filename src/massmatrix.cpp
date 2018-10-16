#include "massmatrix.h"

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  int num_v = F.maxCoeff() + 1;
  Eigen::VectorXd areas(F.rows()), diag(num_v);
  
  igl::doublearea(l, 0., areas);

  for(int i=0; i<F.rows(); i++){
    for(int j=0; j<3; j++){
      // For triangle meshes
      diag(F(i,j)) += areas(i)/6.0;
      // divide by 6 since area is doubled
    }
  }

  M = diag.asDiagonal();
}