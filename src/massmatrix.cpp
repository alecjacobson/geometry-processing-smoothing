#include "massmatrix.h"
#include <math.h>

double triangle_area(double a, double b, double c){
    double s = (a + b + c) / 2.0;
    return sqrt(s * (s-a) * (s-b) * (s-c));
}

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  for (int faceIndex = 0; faceIndex < F.rows(); faceIndex++){
    auto vertices = F.row(faceIndex);
    // Indices
    int v1 = vertices[0];
    int v2 = vertices[1];
    int v3 = vertices[2];

    auto lengths = l.row(faceIndex);
    // Side lengths
    double s1 = lengths[0];
    double s2 = lengths[1];
    double s3 = lengths[2];
    double area = triangle_area(s1, s2, s3);

    // Add mass for each vertex touching the face
    M.diagonal()[v1] += area/ 3.0;
    M.diagonal()[v2] += area/ 3.0;
    M.diagonal()[v3] += area/ 3.0;
  } 
}

