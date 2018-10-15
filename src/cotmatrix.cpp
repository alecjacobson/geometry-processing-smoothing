#include "cotmatrix.h"
#include <Eigen/Core>

double cotangent(double a, double b, double c) {
    double c_half = 1.0 * (a + b + c) / 2;
		double area = sqrt(c_half * (c_half - a) * (c_half - b)* (c_half - c));
    double r_double = 1.0 * (a*b*c) / (2 * area);
    double sin = c / r_double;
	  double cos = 1.0 * (a*a + b*b - c*c) / (2 * a*b);
    double result = cos / (2.0 * sin);
	  return result;
}

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
  // Add your code here
  int V =  F.maxCoeff()+1;
  L.resize(V,V);

  typedef Eigen::Triplet<double> T;
  std::vector<T> tlist;
  tlist.reserve(V*9);

  for (int i=0; i < F.rows(); i++) {
    double l0,l1,l2;
    l2 = l(i,2);
    l0 = l(i,0);
    l1 = l(i,1);
    double ct0,ct1,ct2;
    ct2 = cotangent(l0,l1,l2);
    ct0 = cotangent(l1,l2,l0);
    ct1 = cotangent(l2,l0,l1);
    
    tlist.push_back(T(F(i,0),F(i,1),ct2));
    tlist.push_back(T(F(i,1),F(i,0),ct2));
    tlist.push_back(T(F(i,1),F(i,2),ct0));
    tlist.push_back(T(F(i,2),F(i,1),ct0));
    tlist.push_back(T(F(i,2),F(i,0),ct1));
    tlist.push_back(T(F(i,0),F(i,2),ct1));

    tlist.push_back(T(F(i,0),F(i,0),-(ct2+ct1)));
    tlist.push_back(T(F(i,1),F(i,1),-(ct2+ct0)));
    tlist.push_back(T(F(i,2),F(i,2),-(ct0+ct1)));
  }
  L.setFromTriplets(tlist.begin(),tlist.end());
}

