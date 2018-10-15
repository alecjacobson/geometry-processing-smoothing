#include "cotmatrix.h"

void computeCot(const Eigen::RowVector3d& e, Eigen::RowVector3d& cot) {
  double area = 0.25*std::sqrt(
      (e(0) + e(1) + e(2))*(-e(0) + e(1) + e(2))*
      (e(0) - e(1) + e(2))*(e(0) + e(1) - e(2))
    );

  double sin0 = 2*area/(e(1)*e(2));
  double sin1 = 2*area/(e(0)*e(2));
  double sin2 = 2*area/(e(1)*e(0));

  double cos0 = 0.5*(e(1)*e(1) + e(2)*e(2) - e(0)*e(0))/(e(1)*e(2));
  double cos1 = 0.5*(e(0)*e(0) + e(2)*e(2) - e(1)*e(1))/(e(0)*e(2));
  double cos2 = 0.5*(e(1)*e(1) + e(0)*e(0) - e(2)*e(2))/(e(1)*e(0));

  cot(0) = cos0/sin0;
  cot(1) = cos1/sin1;
  cot(2) = cos2/sin2;
}

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
  // Add your code here
  typedef Eigen::Triplet<double> T;

  // number of vertices
  int v_num = F.maxCoeff() + 1;
  
  // set size of L
  L.resize(v_num, v_num);

  // Build L with triples
  std::vector<T> tripletList;
  tripletList.reserve(F.rows()*3 + v_num);

  // fill off-diagonals
  Eigen::RowVector3d cot;
  for (int i = 0; i < F.rows(); ++i) {
    // compute cotangents
    computeCot(l.row(i), cot);

    // fill L
    const Eigen::RowVector3i& f = F.row(i);

    tripletList.push_back(T(f(0), f(1), 0.5*cot(2)));
    tripletList.push_back(T(f(1), f(2), 0.5*cot(0)));
    tripletList.push_back(T(f(2), f(0), 0.5*cot(1)));

    tripletList.push_back(T(f(1), f(0), 0.5*cot(2)));
    tripletList.push_back(T(f(2), f(1), 0.5*cot(0)));
    tripletList.push_back(T(f(0), f(2), 0.5*cot(1)));

  } // end loop i

  // set from triplets
  L.setFromTriplets(tripletList.begin(), tripletList.end());

  // fill diagonals
  for (int i = 0; i < v_num; ++i) {
    L.coeffRef(i,i) -= L.row(i).sum();
  } // end loop i
  
}

