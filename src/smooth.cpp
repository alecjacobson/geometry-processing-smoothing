#include "smooth.h"
#include "cotmatrix.h"
#include "massmatrix.h"
#include <igl/edge_lengths.h>
#include <iostream>
void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
  // Replace with your code
  Eigen::DiagonalMatrix<double,Eigen::Dynamic>  M;
  Eigen::MatrixXd l;
  Eigen::SparseMatrix<double> L;
  //std::cout << "1" << std::endl;
  igl::edge_lengths(V, F, l);
  cotmatrix(l, F, L);
  //std::cout << "2" << std::endl;
  massmatrix(l, F, M);
  //std::cout << "3" << std::endl;
  Eigen::MatrixXd left = M * G;
  Eigen::MatrixXd temp(M);
  Eigen::SparseMatrix<double> right =  temp.sparseView();
  right = right - lambda * L;
  //std::cout << "4" << std::endl;
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> cg;
  cg.compute(right);
  U = cg.solve(left);

}
