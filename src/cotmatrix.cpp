#include "cotmatrix.h"
#include <iostream>
void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
  // Add your code here
  int num_vertex = F.maxCoeff() + 1;

  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;

  for (int t = 0; t < F.rows(); t++){
    int i = F(t,0);
    int j = F(t,1);
    int k = F(t,2);

    double a = l(t, 0);
    double b = l(t, 1);
    double c = l(t, 2);

    double cos_aij = ((a*a) + (b*b) - (c*c))/(2*a*b);
    double cos_aik = ((c*c) + (b*b) - (a*a))/(2*c*b);
    double cos_ajk = ((c*c) + (a*a) - (b*b))/(2*c*a);

    double sin_aij = sqrt(1 - (cos_aij * cos_aij));
    double sin_aik = sqrt(1 - (cos_aik * cos_aik));
    double sin_ajk = sqrt(1 - (cos_ajk * cos_ajk));

    double cot_aij = cos_aij / sin_aij;
    double cot_aik = cos_aik / sin_aik;
    double cot_ajk = cos_ajk / sin_ajk;

    double A = (.5) * cot_aij;
    double B = (.5) * cot_aik;
    double C = (.5) * cot_ajk;

    tripletList.push_back(T(i,j,A));
    tripletList.push_back(T(j,i,A));
    tripletList.push_back(T(i,k,B));
    tripletList.push_back(T(k,i,B));
    tripletList.push_back(T(j,k,C));
    tripletList.push_back(T(k,j,C));
    tripletList.push_back(T(i,i,-A-B));
    tripletList.push_back(T(j,j,-A-C));
    tripletList.push_back(T(k,k,-B-C));
  }

  L.resize(num_vertex, num_vertex);
  L.setFromTriplets(tripletList.begin(), tripletList.end());

}
