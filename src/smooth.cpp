#include "smooth.h"

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
    // Lengths
    Eigen::MatrixXd l;
    igl::edge_lengths(V, F, l);

    // Cot matrix
    Eigen::SparseMatrix<double> L;
    cotmatrix(l, F, L);

    // Mass matrix
    Eigen::DiagonalMatrix<double,Eigen::Dynamic> M;
    massmatrix(l, F, M);

    // Setup M - lambda*L
    Eigen::SparseMatrix<double> A = -lambda * L;
    for (int i = 0; i < A.rows(); i++){
      A.coeffRef(i, i) += M.diagonal()[i];
    }

    // Cholesky Solver
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> cholesky;
    cholesky.compute(A);
    U = cholesky.solve(M*G);
}