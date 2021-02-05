#include "smooth.h"
#include "massmatrix.h"
#include "cotmatrix.h"
#include "igl/edge_lengths.h"
#include <iostream>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>

typedef Eigen::Triplet<double> T;

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
  Eigen::MatrixXd e;
	Eigen::SparseMatrix<double> L, A;
	Eigen::DiagonalMatrix<double, Eigen::Dynamic> M;

	int v = F.maxCoeff() + 1;

	igl::edge_lengths(V, F, e);
	cotmatrix(e, F, L);
	massmatrix(e, F, M);

	Eigen::MatrixXd temp(v,v);
	temp.setZero();
	temp = M.toDenseMatrix();//Convert to dense matrix for accessing

	A.resize(v, v);
	A = -lambda*L;

  for (int i = 0; i < temp.rows(); i++)
  {
    A.coeffRef(i,i) +=  temp(i,i);
  }

	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
	solver.compute(A);

	U = solver.solve(M * G);

}
