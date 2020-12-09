#include "smooth.h"
#include "massmatrix.h"
#include "cotmatrix.h"
#include <Eigen/SparseCholesky>
#include <igl/edge_lengths.h>

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
    //edge lengths, mass matrix, cotangent matrix
    Eigen::MatrixXd l;
    igl::edge_lengths(V, F, l);
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> M;
    massmatrix(l, F, M);
    Eigen::SparseMatrix<double> L;
    cotmatrix(l, F, L);
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> S;
    
    //A = M - (lambda)*L
    Eigen::SparseMatrix<double> A = -lambda * L;
    for (int i = 0; i < M.rows(); i++) {
        A.coeffRef(i, i) += M.diagonal()[i];
    }
    
    S.compute(A);
    //LDLT solve
    U = S.solve(M * G);
}
