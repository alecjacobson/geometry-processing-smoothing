#include "smooth.h"
#include "cotmatrix.h"
#include "massmatrix.h"
#include <ppl.h>
#include <iostream>
#include <Eigen/Sparse>
#include <concurrent_vector.h>
#include <igl\edge_lengths.h>

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
  // Replace with your code
  int face_num = F.rows();
  int vertex_num = F.maxCoeff() + 1;
  Eigen::MatrixXd l(face_num, 3); l.setZero();
  igl::edge_lengths(V, F, l);

  Eigen::DiagonalMatrix<double, Eigen::Dynamic> M(vertex_num);
  Eigen::SparseMatrix<double>SM(vertex_num, vertex_num); SM.setZero();
  massmatrix(l, F, M);

  auto diag_to_sparse = [&M, &SM](Eigen::DiagonalMatrix<double, Eigen::Dynamic> M,
    Eigen::SparseMatrix<double>&SM)
  {
    typedef Eigen::Triplet<double> tuple;
    Concurrency::concurrent_vector<tuple> tuple_list;
    const int n = M.rows();
    Concurrency::parallel_for(size_t(0), size_t(n), [&M, &SM, &tuple_list](const int m)
    {
      tuple_list.push_back(tuple(m, m, M.diagonal().array()[m]));
    }, Concurrency::auto_partitioner());
    SM.setFromTriplets(tuple_list.begin(), tuple_list.end());
  };

  diag_to_sparse(M, SM);
  Eigen::SparseMatrix<double>L(vertex_num, vertex_num); L.setZero();
  cotmatrix(l, F, L);

  Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver(SM - lambda*L);
  U = solver.solve(M*G);
}