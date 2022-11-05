#include "massmatrix.h"
#include <cmath>

double compute_area_herons_formula(double a, double b, double c);

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  // Add your code here
  int num_vertices = F.maxCoeff() + 1;
  
  Eigen::VectorXd M_vec(num_vertices);
  M_vec.setZero();

  // iterate over the faces
  for (int f = 0; f < F.rows(); f++) {

    // get the lengths
    double a = l(f, 0);
    double b = l(f, 1);
    double c = l(f, 2);

    // get the area of the triangle
    double area = compute_area_herons_formula(a, b, c);

    M_vec(F(f, 0)) += 0.333 * area;
    M_vec(F(f, 1)) += 0.333 * area;
    M_vec(F(f, 2)) += 0.333 * area;
  }

  M.resize(num_vertices);
  M = M_vec.asDiagonal();
}


double compute_area_herons_formula(double a, double b, double c) {

  double semi_per = (a+b+c)/2.0;
  double area_sq = semi_per * (semi_per - a) * (semi_per - b) * (semi_per - c);
  double area = sqrt(area_sq);

  return area;
}

