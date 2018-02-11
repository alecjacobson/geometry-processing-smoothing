#include "massmatrix.h"

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  int num_faces = F.rows();
  int V = F.maxCoeff();
  Eigen::VectorXd masses(V);
  // TODO: can preprocess to do check of is this vertex in triangle faster
  for(int i = 0; i < V; i++){
    // find all triangles t that contain vertex i
    double mass = 0;
    for(int f_idx = 0; f_idx < num_faces; f_idx++){
      if(F(f_idx,0) == i || F(f_idx, 1) == i || F(f_idx, 2) == i){
        // get edge lengths
        double a = l(f_idx, 0); double b = l(f_idx, 1); double c = l(f_idx, 0);
        // Herons formula for area of triangle
        double sp = (a + b + c) / 2; // semiperimeter
        double area = sqrt(sp * (sp - a) * (sp - b) * (sp * c));
        mass += area;
      }
    }
    masses(i) = mass;
  }

  M = masses.asDiagonal();
}
