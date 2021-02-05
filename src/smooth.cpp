#include "smooth.h"
#include <iostream>
#include <igl/edge_lengths.h>
#include <igl/doublearea.h>
#include "cotmatrix.h"
#include "massmatrix.h"
#include<Eigen/SparseCholesky>
#include<Eigen/IterativeLinearSolvers>

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
  // Replace with your code

  // Calculate edge lengths
  Eigen::MatrixXd l;
  igl::edge_lengths(V, F, l);

  // Calculate cotmatrix
  Eigen::SparseMatrix<double> L;
  cotmatrix(l, F, L);
  
  // Calculate massmatrix
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> M;
  massmatrix(l, F, M);

  // Calculate M - lambda * L
  Eigen::SparseMatrix<double> A;
  Eigen::SparseMatrix<double> Md(F.maxCoeff() + 1, F.maxCoeff() + 1);
  Md.setZero();
  std::vector<Eigen::Triplet<double>> tripletsd;
  tripletsd.reserve(F.rows());
  for (int f = 0; f < L.rows(); f++) {
	  tripletsd.push_back({ f, f,  M.diagonal()[f] });
  }
  Md.setFromTriplets(tripletsd.begin(), tripletsd.end());
  A.resize(L.rows(), L.cols());
  A = Md - lambda * L;

  // Solve the linear system
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> sparseSolver;
  sparseSolver.compute(A);
  U = sparseSolver.solve(M * G);



}
