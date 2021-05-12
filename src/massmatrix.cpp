#include "massmatrix.h"
#include <igl/doublearea.h>
void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  // Add your code here
    Eigen::MatrixXd areas;
    igl::doublearea(l, areas);
    M.resize(F.maxCoeff() + 1);
    for (int row = 0; row < F.rows(); ++row){
        for (int col = 0; col < F.cols(); ++col){
            int i = F(row, col);
            M.diagonal()[i] += 1/3.0 * 0.5 * areas(row);
        }
    }

}

