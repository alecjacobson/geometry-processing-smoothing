#include "massmatrix.h"
#include "igl/doublearea.h"

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
    M.resize(F.maxCoeff()+1);
    M.setZero();
    Eigen::MatrixXd area;
    //compute double the area
    igl::doublearea(l, area);
    for (int i = 0; i < F.rows(); i++) {
        //add 1/3 of the area to each
        M.diagonal()(F(i, 0)) += area(i) / 2.0 / 3;
        M.diagonal()(F(i, 1)) += area(i) / 2.0 / 3;
        M.diagonal()(F(i, 2)) += area(i) / 2.0 / 3;
    }
}
