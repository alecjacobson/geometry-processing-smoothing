#include "smooth.h"
#include "igl/edge_lengths.h"
void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
    Eigen::MatrixXd lengths;
    Eigen::SparseMatrix<double> L;
    Eigen::DiagonalMatrix<double,Eigen::Dynamic> M;
    //Need to compute edge-lengths matrix
    
    //Compute Laplacian
    
    //Compute Mass Matrix
    

}
