#include "cotmatrix.h"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <cmath>
// Construct the “cotangent Laplacian” for a mesh with edge lengths l. Each entry
// in the output sparse, symmetric matrix L is given by:
//
// Lij=⎧⎩⎨⎪⎪½cotαij+½cotβij∑j≠iLij0if edge ij existsif i=jotherwise
//
//
//
// Using notation from: https://en.wikipedia.org/wiki/Law_of_cosines
// law of sines and Heron's formula
//
// Triangle with edges a,b,c
//               angles A,B,C
//               indexes 0,1,2 (edge index)
//               angle A is opposite edge a etc.
//
// We aren't allowed to use any trig functions so we will compute
// the cotangent for each triangle angle as follows.
//
//
// cot(A) = cos(A)/sin(A)
// 
// sin(a) = a/d  -- where d is the circumcircle diameter given by:
// d = abc/2A    -- where A is the triangle Area given by Heron's:
// A = sqrt( s(s-a)(s-b)(s-c) ) -- where s is the semiperimeter
// s = (a+b+c)/2
//
// now cos(A) can be found from the law of cosines:
// a^2 = b^2 + c^2 - 2bc cos(A)
// reduces to:
// cos(A) = (b^2 + c^2 -a^2) / 2bc
//
// Now each triangle has 3 angles and 3 contributions across the mesh
// for Lij
// So we will compute them all - assuming a closed mesh with no holes...
// otherwise we'd have to find twinned edges and take them into consideration
//
// sum of each row/col should be 0 -- use that for a check.
// Lij = 1/2 cot( alpha_ij ) + 1/2 cot( beta_ij ) for i != j for each ij edge
// setFromTriplets adds up all the values that have the same coords in the mtx
// convenient...
//


void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
    static bool debug( false );
    
    using Trip =  Eigen::Triplet<double>;
    std::vector< Trip > trips;
    
    const int nTris( F.rows() );
    for( int t=0; t<nTris; ++t )
    {
        Eigen::RowVector3i tri = F.row( t );
        Eigen::RowVector3d e = l.row( t ); // edge lengths - opposite verts
        if( debug )
            std::cout<<"e: "<<e<<std::endl;

        // angle A corresponds is opposite edge 0 (the first edge stored).
         
        // edge lengths opposite verts and angles
        // repeat for all the angles in the triangle
        for( int q=0; q<3; ++q )
        {
            double a = e( q );
            double b = e( (q+1)%3 );
            double c = e( (q+2)%3 );
        
            double s = (a+b+c) / 2;
            double A = std::sqrt( s*(s-a)*(s-b)*(s-c) );
            double d = (a*b*c) / (2*A);

            double sinA = a/d;
            double cosA = ( b*b + c*c - a*a ) / (2*b*c);

            double cotA_half = 0.5*( cosA / sinA );

            // edge coords ij is verts: tri(2)-tri(1) || tri(1)-tri(2)
            int i = tri( (q+1)%3 );
            int j = tri( (q+2)%3 );
            if( debug )
                std::cout<<"t: "<<t<<" q: "<<q<<" i: "<<i<<" j: "<<j<<" val: "
                         <<cotA_half<<std::endl;
            trips.push_back( Trip( i, j, cotA_half ) );
            trips.push_back( Trip( j, i, cotA_half ) );
            // diagonal
            trips.push_back( Trip( i, i, -cotA_half ) );
            trips.push_back( Trip( j, j, -cotA_half ) );
        }
    }

    if( debug )
        std::cout<<"setting triplets: "<<trips.size()<<std::endl;

    const int biggest( F.maxCoeff() );
    L.resize( biggest+1, biggest+1 );  
    L.setFromTriplets( trips.begin(), trips.end() );

    if( debug )
        std::cout<<"L:"<<std::endl<<L<<std::endl;
}
