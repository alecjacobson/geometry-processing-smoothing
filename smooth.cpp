#include "smooth.h"
#include "massmatrix.h"
#include "cotmatrix.h"
#include <igl/edge_lengths.h>
#include <Eigen/SparseCholesky>

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U) {

    // compute the edge lengths
    Eigen::MatrixXd l;
    igl::edge_lengths(V, F, l);
    
    // construct the mass matrix
    Eigen::DiagonalMatrix<double,Eigen::Dynamic> M;
    massmatrix(l, F, M);
      
    // construct the cotangent matrix
    Eigen::SparseMatrix<double> L;
    cotmatrix(l, F, L);
    
    // construct the A matrix in the linear system
    Eigen::SparseMatrix<double> A(M.rows(), M.cols());
    typedef Eigen::Triplet<double> T;
    std::vector<T> diagonal;
    
    Eigen::VectorXd md = M.diagonal();
    for (int i = 0; i < md.size(); i++) {
        diagonal.push_back(T(i,i,md(i)));
    }
    
    A.setFromTriplets(diagonal.begin(), diagonal.end());
    
    A -= lambda*L;
    
    // construct the b matrix in the linear system
    Eigen::MatrixXd b = M*G;
    
    // solve the linear system
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    U = solver.solve(b);
}
