#include "smooth.h"

#include <igl/edge_lengths.h>
#include <cotmatrix.h>
#include <massmatrix.h>

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
    U = G;

    Eigen::MatrixXd l;
    igl::edge_lengths(V, F, l);
    Eigen::SparseMatrix<double> L;
    cotmatrix(l, F, L);
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> M;
    massmatrix(l, F, M);

    int nv = F.maxCoeff() + 1;
    Eigen::MatrixXd b(nv, G.cols());
    Eigen::SparseMatrix<double> A(nv, nv);

    // Fill sparse A with diagonal M
    std::vector<Eigen::Triplet<double>> triplets;
    for (int i = 0; i < M.rows(); ++i) {
        triplets.emplace_back(i, i, M.diagonal()[i]);
    }
    A.setFromTriplets(triplets.begin(), triplets.end());

    b = M * U;
    A = A - lambda*L;

    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
    U = solver.compute(A).solve(b);

}
