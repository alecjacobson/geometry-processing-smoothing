#include "smooth.h"
#include "massmatrix.h"
#include "cotmatrix.h"
#include <igl/edge_lengths.h>

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
    Eigen::MatrixXd l;
    igl::edge_lengths(V, F, l);

    Eigen::SparseMatrix<double> L;
    cotmatrix(l, F, L);

    Eigen::DiagonalMatrix<double, Eigen::Dynamic> M;
    massmatrix(l, F, M);

    Eigen::SparseMatrix<double> A = -lambda * L;
    for (int i = 0; i < M.rows(); ++i) {
        A.coeffRef(i, i) += M.diagonal()(i);
    }

    Eigen::MatrixXd B = M * G;
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver(A);
    U = solver.solve(B);
}
