#include "smooth.h"
#include "igl/edge_lengths.h"
#include "massmatrix.h"
#include "cotmatrix.h"
#include<Eigen/SparseCholesky>
#include "iostream"

void smooth(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        const Eigen::MatrixXd &G,
        double lambda,
        Eigen::MatrixXd &U) {

    //Computing edges lengths
    //columns correspond to edges [1,2],[2,0],[0,1]
    Eigen::MatrixXd l;
    igl::edge_lengths(V, F, l);

    Eigen::SparseMatrix<double> L;
    cotmatrix(l, F, L);

    Eigen::DiagonalMatrix<double, Eigen::Dynamic> M;
    massmatrix(l, F, M);
    // ---

    //Let's create the spatial matrix (M-lambda*L)
    //Unfortunately, Eigen doesnt support this operation (M-lambda*L) in one line. so..
    Eigen::SparseMatrix<double> A = -lambda * L;
    for (int idx = 0; idx < M.rows(); idx++) {
        A.coeffRef(idx, idx) += M.diagonal()(idx);
    }

    // According to the documentation of eigen
    // https://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html
    // this is the recommended Cholesky solver for very sparse and not too large problems

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(A);
    U = solver.solve(M * G);

}
