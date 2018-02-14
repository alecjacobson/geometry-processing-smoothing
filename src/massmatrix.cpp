#include "massmatrix.h"
#include "igl/doublearea.h"

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  int nV = F.maxCoeff() + 1;
  M.resize(nV);

  Eigen::VectorXd dblA;
  igl::doublearea(l, dblA);

  for(int fi = 0; fi < F.rows(); fi++) {
  	for(int vd = 0; vd < 3; vd++) {
  		int vert = F(fi, vd);
  		M.diagonal()[vert] += dblA[fi] / 6.0;
  	}
  }
}

