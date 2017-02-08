#include "cotmatrix.h"
#include <math.h>
#include <iostream>

// Given side lengths a, b, and c, return the
// cosine of angle gamma (opposite from side with length c)
double cosine(double a, double b, double c) {
  return ( pow(a,2) + pow(b,2) - pow(c,2) ) / (double) (2*a*b);
}

// Using trig identity, return the sine of the angle
// defined by the given cosine
double sine(double cosine) {
  return sqrt(1 - pow(cosine,2)); 
}

// Given side lengths a, b and c, return the cotangent
// of angle gamma (opposite from side with length c)
double cotangent(double a, double b, double c) {
  double cos = cosine(a, b, c);
  double sin = sine(cos);
  return cos / sin;
}

void setValueInTripletList(int vi, int vj, double cot, std::vector<Eigen::Triplet<double>> & tripletList) {
  tripletList.push_back(Eigen::Triplet<double>(vi, vj, cot));
  tripletList.push_back(Eigen::Triplet<double>(vj, vi, cot));
  tripletList.push_back(Eigen::Triplet<double>(vi, vi, -1 * cot));
  tripletList.push_back(Eigen::Triplet<double>(vj, vj, -1 * cot));
}

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
  // Get number of vertices
  int v = F.maxCoeff() + 1;
  
  // Construct L
  L.reserve(v*3);
  L.setZero();
  
  // Get all 1/2 cot. Use a tripletList, since this will add up
  // the dupes
  int i, vi, vj;
  double cot;
  std::vector<Eigen::Triplet<double> > tripletList;
  for (i = 0; i < F.rows(); i++) {
      // angle at v0 => opposite of edge v1v2
      vi = F(i,1);
      vj = F(i,2);
      cot = 0.5 * cotangent(l(i,1), l(i,2), l(i,0));
      setValueInTripletList(vi, vj, cot, tripletList);

      // angle at v1 => opposite of edge v2v0
      vi = F(i,2);
      vj = F(i,0);
      cot = 0.5 * cotangent(l(i,2), l(i,0), l(i,1));
      setValueInTripletList(vi, vj, cot, tripletList);
      
      // angle at v2 => opposite of edge v0v1
      vi = F(i,0);
      vj = F(i,1);
      cot = 0.5 * cotangent(l(i,0), l(i,1), l(i,2));
      setValueInTripletList(vi, vj, cot, tripletList);
  }
  L.setFromTriplets(tripletList.begin(), tripletList.end());
}