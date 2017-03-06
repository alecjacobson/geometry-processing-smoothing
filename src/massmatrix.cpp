#include "massmatrix.h"

#include <iostream>

// Construct the diagonal(ized) mass matrix M for a mesh with given face indices
// in F and edge lengths l.
//
// Mij if i = j -->  1/3 Sum_t=1_to_m { Area(t) for triangle t with vertex i
//                                      0 otherwise
//
// so we just have diagnoal values, zero everywhere else
void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
    static bool debug( false );

    const int nTris( F.rows() );
    const int sz( F.maxCoeff() + 1 );

    Eigen::VectorXd myDiag( sz, 1 );
    myDiag.setZero();

    // Would have liked to have used igl::doublearea here but we need V...
    for( int t=0; t<nTris; ++t )
    {
        Eigen::RowVector3i tri = F.row( t );
        Eigen::RowVector3d e = l.row( t ); // edge lengths - opposite verts

        double a = e( 0 );
        double b = e( 1 );
        double c = e( 2 );
        
        double s = (a+b+c) / 2;

        // would be better to do divide by 3 only once along the diag...
        double Ad3 = std::sqrt( s*(s-a)*(s-b)*(s-c) ) / 3; 

        for( int q=0; q<3; ++q )
        {
            int vertId = tri( q );
            myDiag( vertId ) += Ad3;
        }
    }

    M.resize( sz );
    M.setZero();
    M.diagonal() = myDiag;

    if( debug )
        std::cout<<"M is: "<<M.diagonal()<<std::endl;
    
}

