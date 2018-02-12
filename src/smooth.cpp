#include "smooth.h"
#include "igl/edge_lengths.h"
#include "cotmatrix.h"
#include "massmatrix.h"

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
    Eigen::MatrixXd lengths;
    Eigen::SparseMatrix<double> L,A;
    Eigen::DiagonalMatrix<double,Eigen::Dynamic> M;
    Eigen::MatrixXd B;
    
    
    //Need to compute edge-lengths matrix
    lengths.resizeLike(F);
    igl::edge_lengths(V,F,lengths);
    
    
    //Compute Laplacian
    cotmatrix(lengths,F,L);
    //Compute Mass Matrix
    massmatrix(lengths,F,M);
    
    //B matrix
    B.resizeLike(G);
    B = M * G;
    L = lambda * L;
    //A matrix
    A.resize(M.rows(), M.cols());
    A = -L;
    Eigen::VectorXd diagVals = M.diagonal();
    for (int i = 0; i < M.rows(); i ++) {
        A.coeffRef(i,i) += diagVals(i);
    }
    
    U.resizeLike(G);
    
    //Solves for the updates
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
    
    solver.compute(A);
    U = solver.solve(B);
}
