#include "cotmatrix.h"
#include <cmath>

using namespace Eigen;
using namespace std;

double compute_cot(
  double a,
  double b,
  double c)
{
  // compute area: Heron's formula
  double s = (a + b + c) / 2.0;
  double area = sqrt(s * (s - a) * (s - b) * (s - c));
  // cotA = cosA / sinA
  double sinA = 2.0 * area / (b * c);
  double cosA = (b * b + c * c - a * a) * 1.0 / (2 * b * c);
  double cotA = cosA / sinA;
  return cotA;
}

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
  vector<Triplet<double>> triplets;
  for (int i = 0; i < F.rows(); i++) {
    // vertices
    int VA = F(i, 0);
    int VB = F(i, 1);
    int VC = F(i, 2);
    // edges l: columns correspond to edges 23,31,12
    double a = l(i, 0); // BC
    double b = l(i, 1); // AC
    double c = l(i, 2); // CA
    // cot
    double cotA = compute_cot(a, b, c);
    double cotB = compute_cot(b, c, a);
    double cotC = compute_cot(c, a, b);
    // symmetric
    triplets.push_back(Triplet<double>(VA, VB, cotC / 2.0));
    triplets.push_back(Triplet<double>(VB, VA, cotC / 2.0));
    triplets.push_back(Triplet<double>(VB, VC, cotA / 2.0));
    triplets.push_back(Triplet<double>(VC, VB, cotA / 2.0));
    triplets.push_back(Triplet<double>(VC, VA, cotB / 2.0));
    triplets.push_back(Triplet<double>(VA, VC, cotB / 2.0));
    // diagonal
    triplets.push_back(Triplet<double>(VA, VA, -cotC / 2.0));
    triplets.push_back(Triplet<double>(VB, VB, -cotC / 2.0));
    triplets.push_back(Triplet<double>(VB, VB, -cotA / 2.0));
    triplets.push_back(Triplet<double>(VC, VC, -cotA / 2.0));
    triplets.push_back(Triplet<double>(VA, VA, -cotB / 2.0));
    triplets.push_back(Triplet<double>(VC, VC, -cotB / 2.0));
  }
  L.resize(F.maxCoeff() + 1, F.maxCoeff() + 1);
  L.setFromTriplets(triplets.begin(), triplets.end());
}
