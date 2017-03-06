#include "massmatrix.h"

double area(double a, double b, double c)
{
  // calculates area based upon the Heron's formua
  double s = 0.5*(a + b + c);
  return sqrt(s*(s-a)*(s-b)*(s-c));
}

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  /*
    
    construct a mass matrix where each element 
          / ⅓ ∑ⁿₓ₌₁ / Area(x)    if triangle x contains vertex i
    Mᵢⱼ = |         \ 0          otherwise                             if i = j
          \ 0                                                          otherwise

    in simpler terms, the matrix is a diagonal matrix where each entry is the sum of 1/3 of the area of
    all of the triangles attached to the point. 

    Will iterate through all triangles and add 1/3 of that area to the points of triangle
   */

  // creating a diagonal matrix with Eigen we can make a Eigen::VectorXd of #F x 1, and then call .asDiagonal()
  Eigen::VectorXd diagonal = Eigen::VectorXd::Zero(F.maxCoeff() + 1);
  
  for( int32_t i = 0; i < F.rows(); i++)
  {
    double a = area(l(i, 0), l(i, 1), l(i, 2));
    
    diagonal(F(i, 0)) += a/3.0;
    diagonal(F(i, 1)) += a/3.0;
    diagonal(F(i, 2)) += a/3.0;
  }

  M = diagonal.asDiagonal();
}

