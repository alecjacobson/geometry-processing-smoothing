#include "smooth.h"
#include <iostream>
#include <igl/edge_lengths.h>
#include "cotmatrix.h"
#include "massmatrix.h"
#include <Eigen/SparseCholesky>

using namespace std;

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
  int number_of_vertices = V.rows();

  Eigen::MatrixXd edge_lengths;
  igl::edge_lengths(V, F, edge_lengths);

  Eigen::SparseMatrix<double>Laplacian(number_of_vertices, number_of_vertices);
  cotmatrix(edge_lengths, F, Laplacian);

  Eigen::DiagonalMatrix<double,Eigen::Dynamic> Mass(number_of_vertices);
  massmatrix(edge_lengths, F, Mass);

  Eigen::SparseMatrix<double> healthyA;
  // Multiply by the step-size lambda
  healthyA = - lambda * Laplacian;
  // Add in mass, which preserves symmetry
  for (int diagIndex = 0; diagIndex < Mass.rows(); diagIndex++){
    healthyA.coeffRef(diagIndex, diagIndex) += Mass.diagonal()[diagIndex];
  }

  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(healthyA);
  U = solver.solve(Mass * G);

  cout << "Done smoothing!" << endl;
}
