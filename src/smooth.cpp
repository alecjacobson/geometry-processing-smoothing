#include "smooth.h"

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
    Eigen::MatrixXd l;
    igl::edge_lengths(V, F, l);

    Eigen::DiagonalMatrix<double, Eigen::Dynamic> M;
    massmatrix(l, F, M);

    Eigen::SparseMatrix<double> L;
    cotmatrix(l, F, L);

    Eigen::MatrixXd MG = M * G;

    // A = M - lamda * L
    Eigen::SparseMatrix<double> A = -lambda * L;
    for (int i = 0; i < A.rows(); i++) 
    {
        A.coeffRef(i, i) += M.diagonal()[i];
    }

    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> chol(A);
    U = chol.solve(MG);
}
