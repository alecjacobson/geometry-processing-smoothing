#include "cotmatrix.h"
#include <vector>
#include <math.h>

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  
  for (int f = 0; f < F.rows(); f++) {
    for (int v = 0; v < 3; v++) {
      int i = F(f, v);
      int j = F(f, (v+1) % 3);
      
      double eij = l(f, v);
      double ej = l(f, (v+1) % 3);
      double ei = l(f, (v+2) % 3);
      
      double s = (ei + ej + eij) / 2;
      double A = sqrt(s * (s-ei) * (s-ej) * (s-eij));  // heron's formula
      double cot_aij = (ei*ei + ej*ej - eij*eij) / (4 * A);  // law of cosines + sine area formula
      tripletList.push_back(T(i, j, cot_aij / 2));
      tripletList.push_back(T(j, i, cot_aij / 2));
      tripletList.push_back(T(i, i, -cot_aij / 2));
      tripletList.push_back(T(j, j, -cot_aij / 2));
    }
  }
  
  L.setFromTriplets(tripletList.begin(), tripletList.end());
}
