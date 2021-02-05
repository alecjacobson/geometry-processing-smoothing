#include "cotmatrix.h"
#include <iostream>

typedef Eigen::Triplet<double> T;

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
  // setFromTriplets will sum up all entries:
  std::vector<T> triplet_list;
  for (int n = 0; n < F.rows(); n ++) {
    for (int m = 0; m < 3; m ++) {
      int i = F(n, m);
      int j = F(n, (m+1)%3);
      int k = F(n, (m+2)%3);

      double ei = l(n, m);
      double ej = l(n, (m+1)%3);
      double ek = l(n, (m+2)%3);

      double cos_ij = (pow(ei, 2) + pow(ej, 2) - pow(ek, 2)) / (2.0*ei*ej);
      double sin_ij = sqrt(1-pow(cos_ij, 2));
      double cot_ij = 0.5 * (cos_ij / sin_ij);

      triplet_list.push_back(T(i, j, cot_ij));
      triplet_list.push_back(T(j, i, cot_ij));
      triplet_list.push_back(T(i, i, -cot_ij));
      triplet_list.push_back(T(j, j, -cot_ij));
    }
  }
  int num_V = F.maxCoeff() + 1;
  L.resize(num_V, num_V);
  L.setFromTriplets(triplet_list.begin(), triplet_list.end());
}

