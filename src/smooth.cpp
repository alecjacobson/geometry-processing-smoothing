#include "smooth.h"
#include <Eigen/Cholesky>
#include <igl/edge_lengths.h>
#include "massmatrix.h"
#include "cotmatrix.h"

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
  Eigen::MatrixXd l;
  Eigen::SparseMatrix<double> L,A;
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> M;
  igl::edge_lengths(V, F, l);
  cotmatrix(l,F,L);
  massmatrix(l,F,M);
  Eigen::MatrixXd b=M*G;
  typedef Eigen::Triplet<double> T;
  std::vector<T> list;
  Eigen::VectorXd M1=M.diagonal();
  for (int i=0; i<M1.size(); i++)
    list.push_back(T(i,i,M1(i)));
  A.resize(M.rows(),M.rows());
  A.setFromTriplets(list.begin(), list.end()); 
  A-=lambda*L;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(A);
  U=solver.solve(b);
}
