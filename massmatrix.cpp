#include "massmatrix.h"

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M) {
  
    // initialize a vector to store the sums of areas of triangles containing each vertex
    int n = F.maxCoeff() + 1;
    M.resize(n);
    Eigen::VectorXd Area = Eigen::VectorXd::Zero(n);
    
    // loop through each face, compute the area of the triangle at each iteration, and
    // add it to the running total sum of areas of triangles containing each of the relevant
    // vertices
    for (int i = 0; i < F.rows(); i++) {
        
        // get the vertex indices and edge lengths from F and l
        int v0 = F(i, 0);
        int v1 = F(i, 1);
        int v2 = F(i, 2);
        
        double a = l(i,0);
        double b = l(i,1);
        double c = l(i,2);
        
        // Heron's formula for the area of a triangle
        double s = (a + b + c)/2;
        double area = sqrt(s*(s-a)*(s-b)*(s-c));
        
        // add to the running total
        Area(v0) += area;
        Area(v1) += area;
        Area(v2) += area;
    }
    
    // set the diagonal of the mass matrix to (1/3)*A
    M = (1.0/3.0)*Area.asDiagonal();
}

