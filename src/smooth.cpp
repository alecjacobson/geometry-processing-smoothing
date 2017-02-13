#include "smooth.h"
#include "cotmatrix.h"
#include "massmatrix.h"

#include <igl/edge_lengths.h>

#include <Eigen/SparseCholesky>

//
// M^t+1 V^t = (M^t+1 - lambda L^t+1) V^t+1
// assume that small changes in V have a negligeable effect, discretize explicitly:
// M^t V^t = (M^t - lambda L^t) V^t+1
// in our case/notation:
// M^t G^t = (M^t - lambda L^t) U
//
// Solve A x = b
// A -> (M^t - lambda L^t)
// x -> U
// b -> M^t G^t
void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
    static bool debug( false );
    
    Eigen::MatrixXd l;
    igl::edge_lengths( V, F, l );
    
    // get L and M
    Eigen::SparseMatrix<double> L;
    cotmatrix( l, F, L );
    
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> M;
    massmatrix( l, F, M );

    // and our A is (M^t+lambda L^t)
    Eigen::SparseMatrix<double> A = -lambda * L;
    // Have to add M^t to A's diagonal.
    //A = M.transpose() + A
    //A.diagonal() += M.diagonal();
    
    const int n( M.diagonal().size() );
    for( int i=0; i<n; ++i )
    {
        if( debug )
            std::cout<<"A("<<i<<","<<i<<"): "<<A.coeffRef(i,i)<<std::endl;
        A.coeffRef( i, i ) += M.diagonal()(i);
        if( debug )
            std::cout<<"new A("<<i<<","<<i<<"): "<<A.coeffRef(i,i)<<std::endl;
    }

   

    Eigen::MatrixXd b;
    b = M * G;

    if( debug )
        std::cout<<"A:"<<A<<std::endl;
    
    Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > cholesky( A );
    U = cholesky.solve( b );
}
