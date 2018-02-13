#include "massmatrix.h"
#include <igl/doublearea.h>

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  // Add your code here
	
	// compute area for all facets
	Eigen::VectorXd doublearea(F.rows());
	igl::doublearea(l, doublearea);
	Eigen::VectorXd area = 0.5*doublearea;
	// compute diagonal of M as a vector
	int V_count = F.maxCoeff() + 1;
	Eigen::VectorXd diagonal = Eigen::VectorXd::Zero(V_count);
	for (int i = 0; i < F.rows(); i++) {
		diagonal(F(i, 0)) += area(i);
		diagonal(F(i, 1)) += area(i);
		diagonal(F(i, 2)) += area(i);
	}
	diagonal = diagonal / 3.0;
	M.diagonal() = diagonal;
}

